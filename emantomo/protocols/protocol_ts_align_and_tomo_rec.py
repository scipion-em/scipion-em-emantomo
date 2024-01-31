# **************************************************************************
# *
# * Authors:     Scipion Team (scipion@cnb.csic.es)
# *
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import glob
from enum import Enum
from os import remove
from os.path import join, getctime, basename
import numpy as np
import emantomo
from emantomo.constants import TS_DIR, TLT_DIR, TOMOGRAMS_DIR, INTERP_TS, EMAN_ALI_LOSS, ALI_LOSS, TLT_PARAMS, \
    EMAN_OFF_TILT_AXIS
from emantomo.convert import ts2Json, loadJson, convertBetweenHdfAndMrc
from emantomo.objects import EmanMetaData, EmanHdf5Handler
from emantomo.protocols.protocol_base import ProtEmantomoBase, IN_TS
from emantomo.utils import genJsonFileName, getPresentTsIdsInSet
from pwem.convert.headers import fixVolume
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pyworkflow.object import Set, Float
from pyworkflow.protocol import PointerParam, BooleanParam, IntParam, FloatParam, LEVEL_ADVANCED, \
    EnumParam, StringParam, GT, LE
from pwem.protocols import EMProtocol
from pyworkflow.utils import makePath, createLink, Message, replaceExt, getExt
from tomo.objects import SetOfTiltSeries, SetOfTomograms, TiltImage, TiltSeries, Tomogram

# Tomo size choices
SIZE_1K = 0
SIZE_2K = 1
SIZE_4K = 2
R1K = '1k'
R2K = '2k'
R4K = '4k'
OUT_TOMO_SIZE_CHOICES = [R1K, R2K, R4K]
RESOLUTION = {R1K: 1024, R2K: 2048, R4K: 4096}


class outputObjects(Enum):
    tiltSeries = SetOfTiltSeries()
    tiltSeriesInterpolated = SetOfTiltSeries()
    tomograms = SetOfTomograms()


class EmanProtTsAlignTomoRec(ProtEmantomoBase):
    """
    This protocol wraps *e2tomogram.py* EMAN2 program.

    Tomogram tiltseries alignment.
    Tomograms are not normally reconstructed at full resolution, generally limited to 1k x 1k or 2k x 2k,
    but the tilt-series are aligned at full resolution. For high resolution subtomogram averaging, the raw
    tilt-series data is used, based on coordinates from particle picking in the downsampled tomograms.
    On a typical workstation reconstruction takes about 4-5 minutes per tomogram.
    """

    _label = 'TS align & tomo rec'
    _possibleOutputs = outputObjects

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self._tiltAxisAngle = None
        self._finalSamplingRate = None
        self.inTsSet = None
        self._doUpdateTiltAxisAng = False  # If the user introduces a value manually, it must be updated in the metadata

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        alignCond = 'doAlignment'
        recCond = 'doRec'

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TS, PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt Series",
                      important=True,
                      help='Select the set of tilt series to be aligned and/or to reconstruct the corresponding '
                           'tomograms.')
        form.addParam('compressBits', IntParam,
                      default=8,
                      label='Bit compression',
                      expertLevel=LEVEL_ADVANCED,
                      help='Number of bits of precision in outputs with lossless compression. Value -1 means '
                           'uncompressed float')
        form.addParam('keepHdfFile', BooleanParam,
                      default=False,
                      label='Keep hdf files?',
                      expertLevel=LEVEL_ADVANCED,
                      help='If set to Yes, the generated files will be saved in both HDF and MRC formats. They are '
                           'generated in HDF and then converted into MRC. The HDF files are deleted by default to '
                           'save storage.')

        form.addSection(label='TS alignment')
        form.addParam('doAlignment', BooleanParam,
                      default=True,
                      label='Align the tilt series?')
        form.addParam('genInterpolatedTs', BooleanParam,
                      label='Generated interpolated TS?',
                      default=False,
                      condition=alignCond,
                      help='If set to Yes, an additional set of tilt series with the transformation matrix applied '
                           'will be generated. It can be used to check if the alignment was correctly calculated.')
        form.addParam('nLandmarks', IntParam,
                      default=20,
                      condition=alignCond,
                      label='Number of landmarks to use')
        form.addParam('pkKeep', FloatParam,
                      default=0.9,
                      label='Fraction of landmarks to keep in the tracking',
                      validators=[GT(0), LE(1)],
                      expertLevel=LEVEL_ADVANCED,
                      condition=alignCond)
        form.addParam('patchTrack', EnumParam,
                      display=EnumParam.DISPLAY_HLIST,
                      choices=['0', '1', '2'],
                      default=0,
                      condition=alignCond,
                      label='No. patch tracking iterations',
                      expertLevel=LEVEL_ADVANCED,
                      help='use patch tracking before landmark based alignment. input 0/1/2 as the number of patch ' \
                           'tracking iterations.')
        form.addParam('boxSizeTrk', IntParam,
                      default=32,
                      label='Box size of the particles for tracking (pix.)',
                      condition=alignCond,
                      help='It may be helpful to use a larger one for fiducial-less cases.')
        form.addParam('letEmanEstimateTilts', BooleanParam,
                      default=True,
                      label='Should Eman estimate the tilt angles?',
                      condition=alignCond,
                      help='If set to No, a .tlt file will be generated containing the tilt angles read from Scipion '
                           'imported tilt series metadata and passed to Eman. Default=True will let Eman to estimate '
                           'the tilt angles, following the Eman original default behavior.')
        form.addParam('tltStep', FloatParam,
                      default=2,
                      label='Step between tilts',
                      condition=f'{alignCond} and letEmanEstimateTilts')
        form.addParam('zeroId', IntParam,
                      default=-1,
                      label='Index of the center tilt',
                      condition=f'{alignCond} and letEmanEstimateTilts')
        form.addParam('tiltAxisAngle', FloatParam,
                      allowsNull=True,
                      label='Tilt axis angle',
                      condition=alignCond,
                      expertLevel=LEVEL_ADVANCED,
                      help='If not provided, it will be read from the tilt series metadata. If it is not '
                           'contained in the metadata, it will be estimated by EMAN.')
        form.addParam('writeIntermediateRes', BooleanParam,
                      default=False,
                      condition='%s and not genInterpolatedTs' % alignCond,
                      expertLevel=LEVEL_ADVANCED,
                      label='Write intermediate results?',
                      help='They will be generated always the interpolated tilt series are requested.')

        form.addSection(label='Tomogram reconstruction')
        form.addParam('doRec', BooleanParam,
                      default=True,
                      label='Reconstruct the tomograms?')
        form.addParam('outsize', EnumParam,
                      display=EnumParam.DISPLAY_HLIST,
                      default=SIZE_1K,
                      choices=OUT_TOMO_SIZE_CHOICES,
                      condition=recCond,
                      label='Size of output tomograms')
        form.addParam('nIters', StringParam,
                      default='2,1,1,1',
                      condition=recCond,
                      label='No. iterations for bin8, bin4, bin2 images')
        form.addParam('clipz', IntParam,
                      default=-1,
                      condition=recCond,
                      label='Thickness (pix.)',
                      help='Z thickness of the final tomogram output. default is -1, (5/16 of tomogram length).')
        form.addParam('tltkeep', FloatParam,
                      default=0.9,
                      label='Fraction of tilts to keep in the reconstruction',
                      validators=[GT(0), LE(1)])
        form.addParam('bytile', BooleanParam,
                      default=True,
                      condition=recCond,
                      expertLevel=LEVEL_ADVANCED,
                      label='Make final tomogram by tiles?')
        form.addParam('moreTile', BooleanParam,
                      default=False,
                      label='Sample more tiles during rec.?',
                      condition=recCond,
                      expertLevel=LEVEL_ADVANCED,
                      help='If set to Yes, the processing time will be greater, but it can be useful to reduce the '
                           'boundary artifacts when the sample is thick.')
        form.addParam('correctrot', BooleanParam,
                      default=False,
                      label='Correct rotation',
                      condition=recCond,
                      help='Correct for global rotation and position sample flat in tomogram.')
        form.addParam('filterres', FloatParam,
                      default=40,
                      condition=recCond,
                      label='Filter final tomogram to target resolution (Å)')
        form.addParam('rmbeadthr', FloatParam,
                      default=-1.0,
                      label='Density value threshold for removing beads',
                      expertLevel=LEVEL_ADVANCED,
                      condition=recCond,
                      help='"Density value threshold (of sigma) for removing beads. High contrast objects beyond this '
                           'value will be removed. Default is -1 for not removing.')
        form.addParam('extrapad', BooleanParam,
                      default=False,
                      label='Extra pad',
                      expertLevel=LEVEL_ADVANCED,
                      condition=recCond,
                      help='Use extra padding for tilted reconstruction. It is slower and costs more memory, '
                           'but reduces the boundary artifacts when the sample is thick.')
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------

    def _insertAllSteps(self):
        mdObjList = self._initialize()
        for mdObj in mdObjList:
            self._insertFunctionStep(self.convertTsStep, mdObj)
            self._insertFunctionStep(self.writeData2JsonFileStep, mdObj)
            self._insertFunctionStep(self.emanStep, mdObj.tsId)
            self._insertFunctionStep(self.convertOutputStep, mdObj)
            self._insertFunctionStep(self.createOutputStep, mdObj)
        self._insertFunctionStep(self.closeOutputSetsStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.inTsSet = getattr(self, IN_TS).get()
        tltDir = self._getExtraPath(TLT_DIR)
        self.createInitEmanPrjDirs()
        makePath(tltDir)
        if self.tiltAxisAngle.get():  # Value from form
            self._doUpdateTiltAxisAng = True

        mdObjList = []
        presentTsIds = getPresentTsIdsInSet(self.inTsSet)
        tsDict = {ts.getTsId(): ts.clone(ignoreAttrs=[]) for ts in self.inTsSet if ts.getTsId() in presentTsIds}
        counter = 0
        for tsId, ts in tsDict.items():
            if not self.letEmanEstimateTilts.get():
                # Generate tlt file depending on the user choice about the tilt angles
                ts.generateTltFile(join(tltDir, tsId + '.tlt'))
            mdObjList.append(EmanMetaData(tsId=tsId,
                                          ts=ts,
                                          tsHdfName=join(TS_DIR, f'{tsId}.hdf'),
                                          jsonFile=genJsonFileName(self.getInfoDir(), tsId),
                                          processingInd=counter))
            counter += 1
        return mdObjList

    def convertTsStep(self, mObj):
        ts = mObj.ts
        inFile = ts.getFirstItem().getFileName()
        outFile = mObj.tsId + getExt(inFile)
        createLink(inFile, self._getExtraPath(outFile))
        program = emantomo.Plugin.getProgram("e2import.py")
        cmd = [
            f'{outFile}',
            f'--apix={ts.getSamplingRate():.2f}',
            '--import_tiltseries',
            f'--compressbits={self.compressBits.get()}',
            '--importation=copy'
        ]
        self.runJob(program, ' '.join(cmd), cwd=self._getExtraPath())

    def writeData2JsonFileStep(self, mdObj):
        # Generate initial json file required by EMAN for the reconstruction considering a previous alignment
        if self.doRec.get() and not self.doAlignment.get():
            # This only applies to the reconstruction
            ts2Json(mdObj, mode='w')

    def emanStep(self, tsId):
        program = emantomo.Plugin.getProgram("e2tomogram.py")
        cmd = [self.getCommonArgs(tsId)]
        if self.doAlignment.get():
            cmd.append(self._getAlignArgs(tsId))
        if self.doRec.get():
            cmd.append(self._getRecArgs())
        self.runJob(program, ' '.join(cmd), cwd=self._getExtraPath())

    def convertOutputStep(self, mdObj):
        tsId = mdObj.tsId
        hdfTs = self.getCurrentHdfFile(TS_DIR, tsId)
        hdfTomo = self.getCurrentHdfFile(TOMOGRAMS_DIR, tsId)
        filesToConvert = [hdfTs, hdfTomo]
        if self.genInterpolatedTs.get():
            tsHdfInterpFinalLoc = self.getTsInterpFinalLoc(tsId)
            createLink(self._getInterpTsFile(mdObj), tsHdfInterpFinalLoc)
            filesToConvert.append(tsHdfInterpFinalLoc)
        for hdfFile in filesToConvert:
            eh = EmanHdf5Handler(hdfFile)
            extraArgs = '--apix %.3f ' % eh.getSamplingRate()
            convertBetweenHdfAndMrc(self, hdfFile, replaceExt(hdfFile, 'mrc'), extraArgs=extraArgs)
        # Fix the converted tomogram file headerTomograms
        fixVolume(replaceExt(hdfTomo, 'mrc'))

    def createOutputStep(self, mdObj):
        # TS alignment
        if self.doAlignment.get():
            jsonDict = loadJson(mdObj.jsonFile)
            aliLoss = jsonDict[ALI_LOSS]
            alignParams = jsonDict[TLT_PARAMS]
            self._createOutputTs(mdObj, alignParams, aliLoss)
            if self.genInterpolatedTs.get():
                self._createOutputTsInterpolated(mdObj, alignParams)

        # Tomogram reconstruction
        if self.doRec.get():
            tsId = mdObj.tsId
            eh = EmanHdf5Handler(self.getCurrentHdfFile(TOMOGRAMS_DIR, tsId))
            tomogram = Tomogram()
            tomogram.setFileName(self._getOutTomoName(tsId))
            tomogram.setSamplingRate(eh.getSamplingRate())
            tomogram.setTsId(tsId)
            acq = mdObj.ts.getAcquisition()
            if self._doUpdateTiltAxisAng:
                acq.setTiltAxisAngle(self._tiltAxisAngle)
            tomogram.setAcquisition(acq)
            # Set default tomogram origin
            tomogram.setOrigin(newOrigin=False)
            outTomoSet = self.getOutputSetOfTomograms()
            outTomoSet.append(tomogram)
            outTomoSet.update(tomogram)

        self._store()

    def closeOutputSetsStep(self):
        outPutDict = {}
        # TS alignment
        if self.doAlignment.get():
            outTsSet = self.getOutputSetOfTs()
            outTsSet.setStreamState(Set.STREAM_CLOSED)
            outPutDict[outputObjects.tiltSeries.name] = outTsSet
            if self.genInterpolatedTs.get():
                outTsSetInterp = self.getOutputSetOfTs(interpolated=True)
                outTsSetInterp.setStreamState(Set.STREAM_CLOSED)
                outPutDict[outputObjects.tiltSeriesInterpolated.name] = outTsSetInterp

        # Tomogram reconstruction
        if self.doRec.get():
            outTomoSet = self.getOutputSetOfTomograms()
            outTomoSet.setStreamState(Set.STREAM_CLOSED)
            outPutDict[outputObjects.tomograms.name] = outTomoSet

        # Remove HDF files if requested
        # Delete de HDF files if requested
        if not self.keepHdfFile.get():
            hdfExt = '*.hdf'
            tsHdfs = glob.glob(join(self.getTsDir(), hdfExt))
            tomoHdfs = glob.glob(join(self.getTomogramsDir(), hdfExt))
            for hdfFile in tsHdfs + tomoHdfs:
                remove(hdfFile)

        # Define outputs and relations
        self._defineOutputs(**outPutDict)
        if self.doAlignment.get():
            self._defineSourceRelation(self.inputTS.get(), outTsSet)
            if self.genInterpolatedTs.get():
                self._defineSourceRelation(self.inputTS.get(), outTsSetInterp)
        if self.doRec.get():
            self._defineSourceRelation(self.inputTS.get(), outTomoSet)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        errMsg = []
        if not self.doAlignment.get() and not self.doRec.get():
            errMsg.append('No TS alignment nor tomogram reconstruction selected!')
        return errMsg

    # --------------------------- UTILS functions ----------------------------
    def getCommonArgs(self, tsId):
        args = ' %s ' % self._getTsFile(tsId)
        args += '--threads=%i ' % self.numberOfThreads.get()
        return args

        # --------------------------- TS alignment UTILS functions ----------------------------

    def _getTsFile(self, tsId):
        tsFile = glob.glob(self._getExtraPath(TS_DIR, '%s*' % tsId))[0]
        return join(TS_DIR, basename(tsFile))

    def getTAxisAngle(self):
        """The tilt axis angles is expected to be in the metadata at this point, but the user can introduce it manually
        if necessary. The hierarchy to choose the tilt axis angle value is: if no tilt axis value is introduced, it
        will be read from the tilt series metadata."""
        if not self._tiltAxisAngle:
            ts = self.inTsSet.getFirstItem()
            if self.tiltAxisAngle.get():  # Value from form
                self._tiltAxisAngle = self.tiltAxisAngle.get()
            elif ts.getAcquisition().getTiltAxisAngle():  # Value from the TS metadata
                self._tiltAxisAngle = ts.getAcquisition().getTiltAxisAngle()
        return self._tiltAxisAngle

    def _getAlignArgs(self, tsId):
        tAx = self.getTAxisAngle()
        args = [
            f'--npk={self.nLandmarks.get()}',
            f'--pkkeep={self.pkKeep.get():.2f}',
            f'--bxsz={self.boxSizeTrk.get()}',
            f'--patchtrack={self.patchTrack.get()}',
            f'--compressbits={self.compressBits.get()}'
        ]
        if self.letEmanEstimateTilts.get():
            # Introduce angular step and zero id
            args.append(f'--tltstep={self.tltStep.get():.2f}')
            args.append(f'--zeroid={self.zeroId.get()}')
        else:
            args.append(f'--rawtlt={join(TLT_DIR, tsId + ".tlt")}')
        if not self.writeIntermediateRes.get() and not self.genInterpolatedTs.get():
            args.append('--notmp')
        if tAx:
            args.append(f'--tltax={tAx:.2f}')
        args.append('--verbose=9')
        return ' '.join(args)

    def getOutputSetOfTs(self, interpolated=False):
        inTs = self.inputTS.get()
        outTs = self.getTiltSeries(interpolated=interpolated)
        if outTs:
            outTs.enableAppend()
            tiltSeries = outTs
        else:
            if interpolated:
                suffix = 'interpolated'
                x, y, z, _ = ImageHandler().getDimensions(glob.glob(self._getExtraPath(TS_DIR, '*_interp.mrc'))[0])
                dims = (x, y, z)
            else:
                suffix = ''
                dims = inTs.getDim()

            tiltSeries = SetOfTiltSeries.create(self._getPath(), template='tiltseries', suffix=suffix)
            tiltSeries.copyInfo(self.inputTS.get())
            tiltSeries.setDim(dims)
            tiltSeries.setSamplingRate(inTs.getSamplingRate())
            tiltSeries.setStreamState(Set.STREAM_OPEN)
            self.setTiltSeries(tiltSeries, interpolated=interpolated)

        return tiltSeries

    def getTiltSeries(self, interpolated=False):
        if interpolated:
            return self.getObjByName(outputObjects.tiltSeriesInterpolated.name)
        else:
            return self.getObjByName(outputObjects.tiltSeries.name)

    def setTiltSeries(self, tsSet, interpolated=False):
        if interpolated:
            setattr(self, outputObjects.tiltSeriesInterpolated.name, tsSet)
        else:
            setattr(self, outputObjects.tiltSeries.name, tsSet)

    def _createOutputTs(self, mdObj, alignParams, aliLoss):
        outTsSet = self.getOutputSetOfTs()
        tiltSeries = self._createCurrentOutTs(mdObj.ts)
        outTsSet.append(tiltSeries)
        for idx, ti in enumerate(mdObj.ts.iterItems()):
            outTi = TiltImage()
            outTi.copyInfo(ti, copyId=True)
            outTi.setIndex(ti.getIndex())
            outTi.setFileName(self._getExtraPath(TS_DIR, mdObj.tsId + '.mrc'))
            curerntAlignParams = alignParams[idx]
            setattr(outTi, EMAN_ALI_LOSS, Float(aliLoss[idx]))
            self._genTrMatrixFromEmanAlign(outTi, curerntAlignParams)
            # Append the current tilt image to the corresponding tilt series
            tiltSeries.append(outTi)

        outTsSet.update(tiltSeries)
        return outTsSet

    def _createOutputTsInterpolated(self, mdObj, alignParams):
        outTsSet = self.getOutputSetOfTs(interpolated=True)
        tiltSeries = self._createCurrentOutTs(mdObj.ts, interpolated=True)
        outTsSet.append(tiltSeries)
        for idx, ti in enumerate(mdObj.ts.iterItems()):
            outTi = TiltImage()
            finalName = self.getTsInterpFinalLoc(mdObj.tsId).replace('.hdf', '.mrc')
            outTi.copyInfo(ti, copyId=True)
            outTi.setTsId(outTi.getTsId())
            outTi.setIndex(ti.getIndex())
            outTi.setFileName(finalName)
            setattr(outTi, EMAN_OFF_TILT_AXIS, Float(alignParams[idx][4]))  # Off tilt axis angle, extended parameter
            outTi.setTiltAngle(alignParams[idx][3])  # Refined tilt angle
            tiltSeries.append(outTi)

        outTsSet.setAlignment2D()
        outTsSet.update(tiltSeries)
        return outTsSet

    @staticmethod
    def genUpdateTsInterpOrigin(dims, samplingRate):
        origin = Transform()
        x = dims[0]
        y = dims[1]
        origin.setShifts(-x/2 * samplingRate,
                         -y/2 * samplingRate,
                         0)
        return origin

    def _createCurrentOutTs(self, ts, interpolated=False):
        tiltSeries = TiltSeries()
        if interpolated:
            # Avoid the copyInfo and set the origin manually to get the squared TS dims and origin correctly stored
            tsId = ts.getTsId()
            sr = ts.getSamplingRate()
            tiltSeries.setSamplingRate(sr)
            tiltSeries.setAcquisition(ts.getAcquisition())
            tiltSeries.setTsId(tsId)
            x, y, z, _ = ImageHandler().getDimensions(self._getExtraPath(TS_DIR, tsId + '.mrc'))
            dims = (x, y, z)
            tiltSeries.setDim(dims)
            tiltSeries.setOrigin(self.genUpdateTsInterpOrigin(dims, sr))
        else:
            tiltSeries.copyInfo(ts)
        if self._doUpdateTiltAxisAng:
            acq = ts.getAcquisition()
            acq.setTiltAxisAngle(self._tiltAxisAngle)
            tiltSeries.setAcquisition(acq)
        return tiltSeries

    def _getInterpTsFile(self, mdObj):
        matchingFiles = glob.glob(self._getExtraPath('tomorecon_%02d' % mdObj.processingInd, INTERP_TS + '.hdf'))
        return max(matchingFiles, key=getctime)

        # --------------------------- reconstruct tomograms UTILS functions ----------------------------

    def _getRecArgs(self):
        args = [f'--outsize={OUT_TOMO_SIZE_CHOICES[self.outsize.get()]}',
                f'--niter={self.nIters.get()}',
                f'--clipz={self.clipz.get()}',
                f'--filterres={self.filterres.get():.2f}',
                f'--rmbeadthr={self.rmbeadthr.get():.2f}',
                f'--tltkeep={self.tltkeep.get()}']

        if self.doRec.get() and not self.doAlignment.get():
            args.append('--load')  # Load existing tilt parameters
            args.append('--noali')  # Skip initial alignment

        if self.bytile.get():
            args.append('--bytile')
            if self.doAutoclipXY():
                args.append('--autoclipxy')

        if self.moreTile.get():
            args.append('--moretile')

        if self.correctrot.get():
            args.append('--correctrot')

        if self.extrapad.get():
            args.append('--extrapad')

        return ' '.join(args)

    # def getFinalSampligRate(self):
    #     # It has to be recalculated due to EMAN size management
    #     if not self._finalSamplingRate:
    #         tsSet = self.inputTS.get()
    #         sizeX = tsSet.getDimensions()[0]
    #         exponent = np.ceil(np.log2(sizeX / RESOLUTION[OUT_TOMO_SIZE_CHOICES[self.outsize.get()]]).clip(min=0))
    #         binning = 2 ** exponent
    #         self._finalSamplingRate = binning * tsSet.getSamplingRate()
    #     return self._finalSamplingRate

    def _getOutTomoName(self, tsId):
        return glob.glob(self._getExtraPath(TOMOGRAMS_DIR, '%s*.mrc' % tsId))[0]

    def getOutputSetOfTomograms(self):
        outTomograms = self.getTomograms()
        if outTomograms:
            outTomograms.enableAppend()
            tomograms = outTomograms
        else:
            tomograms = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
            tomograms.copyInfo(self.inputTS.get())
            # Get the sampling rate from the header of a HDF file present in the tomograms directory
            eh = EmanHdf5Handler(glob.glob(join(self.getTomogramsDir(), '*.hdf'))[0])
            tomograms.setSamplingRate(eh.getSamplingRate())
            tomograms.setStreamState(Set.STREAM_OPEN)
            self.setTomograms(tomograms)

        return tomograms

    def getTomograms(self):
        return getattr(self, outputObjects.tomograms.name, None)

    def setTomograms(self, tomoSet):
        setattr(self, outputObjects.tomograms.name, tomoSet)

    def getCurrentHdfFile(self, dirName, tsId):
        return glob.glob(self._getExtraPath(dirName, tsId + '*.hdf'))[0]

    def getTsInterpFinalLoc(self, tsId):
        return self._getExtraPath(TS_DIR, tsId + '_interp.hdf')

    @staticmethod
    def _genTrMatrixFromEmanAlign(outTi, alignParams):
        transform = Transform()
        matrix = np.eye(3)
        rotMatrix = np.eye(2)
        # tlt_params: (N, 5) list, where N is the number of tilt in the tilt series. The columns are translation
        # along x,y axis, and rotation around z, y, x axis in the EMAN coordinates. The translation is in unbinned
        # pixels, and rotation is in degrees
        rotAngle = -np.deg2rad(alignParams[2])
        tiltAngleRefined = alignParams[3]
        offTiltAngle = alignParams[4]
        tx = -alignParams[0]
        ty = -alignParams[1]
        rotAngleCorrected = rotAngle + np.deg2rad(offTiltAngle)
        rotMatrix[0, 0] = rotMatrix[1, 1] = np.cos(rotAngleCorrected)
        rotMatrix[0, 1] = np.sin(rotAngleCorrected)
        rotMatrix[1, 0] = -np.sin(rotAngleCorrected)
        shiftsEman = np.array([tx, ty])
        shiftsImod = rotMatrix.dot(shiftsEman)
        matrix[0, 0] = matrix[1, 1] = rotMatrix[0, 0]
        matrix[0, 1] = rotMatrix[0, 1]
        matrix[1, 0] = rotMatrix[1, 0]
        matrix[0, 2] = shiftsImod[0]
        matrix[1, 2] = shiftsImod[1]
        transform.setMatrix(matrix)
        outTi.setTransform(transform)
        outTi.setTiltAngle(tiltAngleRefined)
        setattr(outTi, EMAN_OFF_TILT_AXIS, Float(offTiltAngle))  # Extended parameter

    def doAutoclipXY(self):
        """Code behavior expected for EMAN's option autoclipxy to be automatic.
        From EMAN tutorial https://blake.bcm.edu/emanwiki/EMAN2/e2tomo_p22:
        'Always check --autoclipxy for non-square micrographs'."""
        ih = ImageHandler()
        x, y, _, _ = ih.getDimensions(self.inTsSet.getFirstItem().getFirstItem().getFileName())
        return True if x != y else False
