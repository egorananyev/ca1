####!/usr/bin/arch -i386 /usr/bin/python # -*- coding: utf-8 -*-
"""
Consciousness, attention, and depth of suppression.
2018-02-05
Egor Ananyev
"""

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import visual, core, data, event, gui
#from psychopy import gui
from psychopy.constants import *  # things like STARTED, FINISHED
import numpy as np
import pandas as pd
from datetime import datetime
import os, shutil, itertools, copy  # handy system and path functions
import pyglet
import glob

# Initiating the keyboard
from psychopy.iohub import launchHubServer
io = launchHubServer()
kb_device = io.devices.keyboard

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))

# ====================================================================================
## Initial variables.
projName = 'ca1'
# Window circles (specified in degrees of visual angles [dva]):
#winSz = 7.2 # 5.03; calculated as 5/x=sqrt(2)/2 => x=10/sqrt(2)
winOffX = 4.25 # 6 # 5.62
winOffY = 2.5 # 5.5 (3.5cm ~= 124px)
winThickness = 2 # in pixels
# Timing variables:
ISIduration = .25
fixSz = .15
# MCs:
precompileMode = 1 # get the precompiled MCs
grtSize = 256 # size of 256 is 71mm, or 7.2dova
#contrSteps = [1.3,1.3,.9,.9,.6,.6,.4,.4,.3,.3] #10
#contrSteps = [.2,.2,.1,.1,.05,.05,.02,.02,.01,.01,.005,.005] #10
contrSteps = [.5,.5,.3,.3,.2,.2,.15,.15,.12,.12,.1,.1]
# Dimensions:
###### 7.2dova = 71mm = 256px; 475x296mm, 563mm viewing dist ######
dr = (1680,1050) # display resolution in px
#dr = (1366,768)
#dd = (47.5,29.6) # display dimensions in cm
dd = (29.5,16.6)
ds = 50+2.5+3.5 #49.5 # distance to screen in cm
trialNfb = False # do we give the trial number feedback?

# ====================================================================================
# Converter functions:
def cm2px(cm,dr=dr,dd=dd):
    px = int(cm*(dr[0]/dd[0]))
    return px
def px2cm(px,dr=dr,dd=dd):
    cm = px/(dr[0]/dd[0])
    return cm
def cm2dg(cm,ds=ds):
    dg = np.degrees(np.arctan(cm/ds))
    return dg
def dg2cm(dg,ds=ds):
    cm = ds*np.tan(np.radians(dg))
    return cm
def px2dg(px,cm2dg=cm2dg,px2cm=px2cm):
    dg = cm2dg(px2cm(px))
    return dg
def dg2px(dg,cm2px=cm2px,dg2cm=dg2cm):
    px = int(cm2px(dg2cm(dg)))
    return px

# ====================================================================================
# Converting win dimensions to pixels
#winSz = dg2px(winSz)
winSz = grtSize + 2
winOffX = dg2px(winOffX)
winOffY = dg2px(winOffY)
fixSz = dg2px(fixSz)
posCentL = [-winOffX, winOffY]
posCentR = [winOffX, winOffY]
#print winSz 
#print posCentL 
#print posCentR 

# ====================================================================================
# Store info about the experiment session
expInfo = {u'exp': '', u'subj': u'',  u'phase': '', u'sess': u''}
# experiments:
### 1 = mask contrast manipulation
### 2 = mask speed
### 3 = mask-target orientation match
# dom 0 = left, 1 = right, '' = unkown (do the domTest)
dlg = gui.DlgFromDict(dictionary=expInfo, title=projName,
                      order=['exp','subj','phase','sess']) # dialogue box
if dlg.OK == False: core.quit()  # user pressed cancel
timeNow = datetime.now()
expInfo['time'] = datetime.now().strftime('%Y-%m-%d_%H%M')
# do the dominance test if the paradigm field is left blank or at 'bc3-':
if expInfo['exp']=='' and expInfo['phase']=='': # phase0: domtest
    domTest = True
    expNum = 'dom'
    expPhase = 0
    expInfo['dom'] = -1
else: # phases 1 or 2:
    expPhase = int(expInfo['phase'])
    domTest = False
    expNum = expInfo['exp']
    # the dominant eye is figured out later (as we need some functions later in code)

# ====================================================================================

# Setup the Window
win = visual.Window(size=dr, fullscr=False, screen=1, allowGUI=False, 
      allowStencil=False, color='grey', blendMode='avg', useFBO=True, units='pix')
# store frame rate of monitor if we can measure it successfully:
frameRate=win.getActualFrameRate()
if frameRate!=None:
    frameDur = 1.0/round(frameRate)
else:
    frameDur = 1.0/60.0 # couldn't get a reliable measure so guess

# ====================================================================================

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
if precompileMode:
    precompiledDir = '../../../mc/precompiledMCs'
dataDir = '..' + os.sep + 'data'
if expNum == 'dom':
    fileName = '%s_%s_subj%s_phase%s_sess%s_%s' %(projName, expNum, expInfo['subj'], 
        expInfo['phase'], expInfo['sess'], expInfo['time'])
else:
    fileName = '%s_exp%s_subj%s_phase%s_sess%s_%s' %(projName, expNum, expInfo['subj'], 
        expInfo['phase'], expInfo['sess'], expInfo['time'])
filePath = dataDir + os.sep + fileName
print filePath

# Condition-related variables
if expNum == 'dom':
    condFilePath = 'cond-files'+os.sep+'cond-'+projName+'_'+expNum+'.csv'
else:
    condFilePath = 'cond-files'+os.sep+'cond-'+projName+'_exp'+expNum+'.csv'
print condFilePath
os.chdir(_thisDir)

# ====================================================================================

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Initialize components for Routine "instructions"
instructionsClock = core.Clock()
instrTextL = visual.TextStim(win, text='press any key\n     to start', font='Cambria',
                             pos=posCentL, height=dg2px(.5), wrapWidth=dg2px(4),
                             color='white', alignHoriz='center')
instrTextR = visual.TextStim(win, text='press any key\n     to start', font='Cambria',
                             pos=posCentR, height=dg2px(.5), wrapWidth=dg2px(4),
                             color='white', alignHoriz='center')

# Initialize components for Routine "trial"
trialClock = core.Clock()
moveClock = core.Clock()
maskMoveClock = core.Clock()
ISI = core.StaticPeriod(win=win, screenHz=frameRate, name='ISI')
# circular windows:
winL = visual.Polygon(win, edges=36, size=[winSz, winSz], pos=posCentL,
                      lineWidth=winThickness, lineColor='white')
winR = visual.Polygon(win, edges=36, size=[winSz, winSz], pos=posCentR,
                      lineWidth=winThickness, lineColor='white')
# target gabor (serves as prime in phase 2)
targGab = visual.GratingStim(win, tex='sin', mask='circle', size=[winSz, winSz],
                      pos=posCentL) # this position will change dep. on eyeDom
# target blob
targBlob = visual.Polygon(win, radius=winSz, edges=72, pos=posCentL)
# fixation:
fixL = visual.ShapeStim(win, pos=posCentL, vertices=((0,-fixSz), (0,fixSz), (0,0), 
                                                     (-fixSz,0), (fixSz,0)),
                        lineWidth=.2, closeShape=False, lineColor='white')
fixR = visual.ShapeStim(win, pos=posCentR, vertices=((0,-fixSz), (0,fixSz), (0,0), 
                                                     (-fixSz,0), (fixSz,0)),
                        lineWidth=.2, closeShape=False, lineColor='white')
# Trial number feedback:
if trialNfb:
    trialNfbText = visual.TextStim(win=win, text='', font='Cambria', 
                                   pos=(0,0), height=dg2px(.55), wrapWidth=dg2px(4.5),
                                   color='white')
# pause text:
#vertShift = np.array([0,40])
pauseStr = 'press spacebar\n\n    to continue'
pauseTextL = visual.TextStim(win, text=pauseStr, font='Cambria',
                             alignHoriz='center', pos=(posCentL), height=dg2px(.5),
                             wrapWidth=None, color='white') #dg2px(3.5)
pauseTextR = visual.TextStim(win, text=pauseStr, font='Cambria',
                             alignHoriz='center', pos=(posCentR), height=dg2px(.5),
                             wrapWidth=dg2px(3.5), color='white')

# resp question & feedback:
fdbStrQn = '?'
fdbStrL = '<'
fdbStrR = '>'
fdbStr0 = 'o'
fdbTextL = visual.TextStim(win, text=fdbStrQn, font='Cambria',
                          alignHoriz='center', pos=(posCentL), height=dg2px(.4),
                          wrapWidth=dg2px(3), color='white')
fdbTextR = visual.TextStim(win, text=fdbStrQn, font='Cambria',
                          alignHoriz='center', pos=(posCentR), height=dg2px(.4),
                          wrapWidth=dg2px(3), color='white')

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

#------Prepare to start Routine "instructions"-------
t = 0
instructionsClock.reset()  # clock 
frameN = -1
# update component parameters for each repeat
instrKey = event.BuilderKeyResponse()  # create an object of type KeyResponse
instrKey.status = NOT_STARTED
# keep track of which components have finished
instructionsComponents = []
instructionsComponents.append(instrTextL)
instructionsComponents.append(instrTextR)
instructionsComponents.append(instrKey)
for thisComponent in instructionsComponents:
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED


# ====================================================================================
# Setting up the conditions:
condList = data.importConditions(condFilePath)
stairs, completedStairs = [], []
for thisCond in condList:
    thisInfo = copy.copy(thisCond)
    stairLabel = 'st' + str(thisCond['startContr'])
    if thisCond['targXoff1']<1: targLoc = 'L'
    else: targLoc = 'R'
    if expNum=='dom': stairLabel += '_targEyeR' + str(thisCond['targEyeR'])
    if expNum=='1': stairLabel += '_targLoc' + targLoc + \
                                  '_mcContr' + str(thisCond['mcContr'])
    if expNum=='2': stairLabel += '_targLoc' + targLoc + \
                                  '_mcBv' + str(thisCond['mcBv'])
    if expNum=='3': stairLabel += '_targLoc' + targLoc + \
                                  '_mcBtheta' + str(thisCond['mcBtheta'])
    thisInfo['label'] = stairLabel
    thisStair = data.StairHandler(startVal = thisInfo['startContr'],
                                  extraInfo = thisInfo, maxVal=0, minVal=-2,
                                  nReversals = thisInfo['nRevs'],
                                  nUp = 1, nDown = 1, stepType='lin',
                                  stepSizes = contrSteps[0:thisInfo['nRevs']],
                                  name = stairLabel)
    stairs.append(thisStair)

# An empty data set for storing behavioural responses:
behResp = []
    
# Creating a copy of the Conditions file for book-keeping and analyses:
if not os.path.exists(filePath):
    os.makedirs(filePath)
shutil.copyfile(condFilePath, filePath + os.sep + 
                os.path.basename(condFilePath))
outFileName = filePath + os.sep + fileName + '.csv'
threshFileName = filePath + os.sep + fileName + '_thresh.csv'

# ====================================================================================
# Various functions for use in trials:

def writeStair(thisStair, filePath):
    stairFileName = filePath + os.sep + thisStair.name
    thisStair.saveAsPickle(stairFileName)
    thisStair.saveAsText(stairFileName)

def dfStair(thisStair, projName, expNum, targEyeR):
    mcPeriFade = thisStair.extraInfo['mcPeriFade']
    mcPeriGap = np.int(thisStair.extraInfo['mcSz']/2-mcPeriFade)
    meanRev6 = np.average(thisStair.reversalIntensities[-6:])
    # have the information recorded in a csv file as well:
    dT = pd.DataFrame({'projName': projName, 'expNum': expNum,
                       'time': expInfo['time'],
                       'subj': expInfo['subj'],
                       'dom': expInfo['dom'],
                       'sess': expInfo['sess'],
                       'nRevs': thisStair.extraInfo['nRevs'],
                       'mcSz': thisStair.extraInfo['mcSz'],
                       'mcSf': thisStair.extraInfo['mcSf'],
                       'mcBv': thisStair.extraInfo['mcBv'],
                       'mcBsf': thisStair.extraInfo['mcBsf'],
                       'mcPeriGap': mcPeriGap,
                       'mcPeriFade': mcPeriFade,
                       'targSz': thisStair.extraInfo['targSz'],
                       'targSf': thisStair.extraInfo['targSf'],
                       'targOri1': thisStair.extraInfo['targOri1'],
                       'targOri2': thisStair.extraInfo['targOri2'],
                       'targXoff1': thisStair.extraInfo['targXoff1'],
                       'targXoff2': thisStair.extraInfo['targXoff2'],
                       'targYoff': thisStair.extraInfo['targYoff'],
                       'targV': thisStair.extraInfo['targV'],
                       'targTtot': thisStair.extraInfo['targTtot'],
                       'targTpeak': thisStair.extraInfo['targTpeak'],
                       'targEyeR': targEyeR,
                       'trialT': thisStair.extraInfo['trialT'],
                       'fixCross': thisStair.extraInfo['fixCross'],
                       'stairLabel': thisStair.extraInfo['label'],
                       'stairStart': thisStair.extraInfo['startContr'],
                       'meanRev6': [meanRev6]})
    return dT

dataCols = ['projName', 'expNum', 'time', 'subj', 'dom', 'sess', 'nRevs',
            'mcSz', 'mcSf', 'mcBv', 'mcBsf', 'mcPeriGap', 'mcPeriFade', 
            'targSz', 'targSf', 'targOri1', 'targOri2', 'targXoff1', 'targXoff2',
            'targYoff', 'targV', 'targTtot', 'targTpeak', 'trialT',
            'fixCross', 'stairLabel', 'stairStart', 'meanRev6']

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def sigmoidMod(x): # modified such that 0->0, 1->1
    return 1 / (1 + np.exp(-x*10+5))

x = np.arange(-grtSize/2,grtSize/2)
y = np.arange(-grtSize/2,grtSize/2)
x, y = np.meshgrid(x, y)
R = np.sqrt((x+.5)**2 + (y+.5)**2) # adding .5 ensures symmetry

def periMask(periGap, periFade, R=R):
    return sigmoid(R * (-10./(periFade)) + 5 + periGap*(10./periFade))*2 - 1

# ====================================================================================

## Figure out the dominance for phases 1-2:

def dfThreshDomRaw(dfThresh): # returns threshold averages from dominance test:
    dfThresh['eye'] = [x[-1] for x in dfThresh.index.values.tolist()]
    return dfThresh.groupby('eye')[-6:].mean().mean(axis=1)

def dfThreshDom6(dfThresh): # returns threshold averages from dominance test:
    dfThresh['eye'] = [x[-1] for x in dfThresh['stairLabel'].tolist()]
    groupByEye = dfThresh.groupby('eye')['meanRev6'].mean()
    return groupByEye

def dfThreshPhase1(dfThresh): # returns threshold averages from dominance test:
    dfThresh['loc'] = [x[-12] for x in dfThresh['stairLabel'].tolist()]
    groupByLoc = dfThresh.groupby('loc')['meanRev6'].mean()
    return groupByLoc

# For both 1st and 2nd phase, decide on the target eye:
if expPhase > 0:
    # Take last dom test for subj for data on thresholds:
    domDir = glob.glob(dataDir + os.sep + '*dom_subj' + expInfo['subj'] + '*')[-1]
    domCsv = glob.glob(domDir + os.sep + '*dom*[0-9].csv')[0]
    dfDomThresh = dfThreshDom6(pd.read_csv(domCsv, index_col=False))
    print dfDomThresh
    if dfDomThresh[0] < dfDomThresh[1]: expInfo['dom'] = 0 # left eye dom
    else: expInfo['dom'] = 1 # right eye dom; default if equal thresholds
    targEyeR = 2-(2**int(float(expInfo['dom'])))
    print 'targEye = ' + str(targEyeR)

if expPhase == 2:
    # Take threshold measurements from phase 1
    phase1dir = glob.glob(dataDir + os.sep + '*exp' + expNum + \
                          '_subj' + expInfo['subj'] + '_phase1*')[-1]
    phase1csv = glob.glob(phase1dir + os.sep + '*phase1*[0-9].csv')[0]
    # Assuming that the thresholds for dimmer masks are lower than for brighter masks 
    # (etc for exp2-3), take the averages for each location for mcContr==0.5. The result
    # is 2 averaged thresholds, one for each location.
    dfPhase1thresh = dfThreshPhase1(pd.read_csv(phase1csv, index_col=False))
    print dfPhase1thresh

# ====================================================================================

#-------Start Routine "instructions"-------
continueRoutine = True
while continueRoutine:
    # get current time
    t = instructionsClock.getTime()
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *instrText* updates
    if t >= 0.0 and instrTextL.status == NOT_STARTED:
        # keep track of start time/frame for later
        instrTextL.tStart = t  # underestimates by a little under one frame
        instrTextL.frameNStart = frameN  # exact frame index
        instrTextL.setAutoDraw(True)
        instrTextR.tStart = t  # underestimates by a little under one frame
        instrTextR.frameNStart = frameN  # exact frame index
        instrTextR.setAutoDraw(True)
    
    # *instrKey* updates
    if t >= 0.0 and instrKey.status == NOT_STARTED:
        # keep track of start time/frame for later
        instrKey.tStart = t  # underestimates by a little under one frame
        instrKey.frameNStart = frameN  # exact frame index
        instrKey.status = STARTED
        # keyboard checking is just starting
        event.clearEvents(eventType='keyboard')
        winL.setAutoDraw(True)
        winR.setAutoDraw(True)
    if instrKey.status == STARTED:
        theseKeys = event.getKeys()
        
        # check for quit:
        if "escape" in theseKeys:
            endExpNow = True
        if len(theseKeys) > 0:  # at least one key was pressed
            # a response ends the routine
            continueRoutine = False
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        routineTimer.reset()  # if we abort early the non-slip timer needs reset
        break
    continueRoutine = False  # reverts to True if at least 1 component still running
    for thisComponent in instructionsComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # check for quit (the Esc key)
    if endExpNow or event.getKeys(keyList=["escape"]):
        core.quit()
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()
    else:  # this Routine was not non-slip safe so reset non-slip timer
        routineTimer.reset()

#-------Ending Routine "instructions"-------
for thisComponent in instructionsComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)

# ====================================================================================
# Initiating the stair loop

nTrialsDone = 0
nStairsDone = 0
while len(stairs)>0:

    # printing current reversals:
    curRevs = []
    nRevs = thisStair.extraInfo['nRevs']
    nStairsTotal = len(stairs) + nStairsDone
    for thisStair in stairs:
        curRevs.append(len(thisStair.reversalIntensities))
    #if nStairsDone>0:
    #    curRevs = [curRevs, np.repeat(nRevs,nStairsDone)]
    #    #curRevs.append(np.repeat(10,nStairsDone))
    #curPerc = np.average(curRevs)/np.repeat(nRevs,nStairsTotal)*100
    #print 'staircase progress: %i%%' %(curPerc.item(0))
    print 'curRevs: ' + str(curRevs)

    # selecting a stair for this trial:
    np.random.shuffle(stairs)
    thisStair = stairs.pop()

    try:
        # current contrast:
        thisContr = thisStair.next() # contrast value
        contrStr = 'start=%.1f, cur=%.2f' %(thisStair.extraInfo['startContr'], thisContr)

    except StopIteration:
        print '-------------------------------------------------'
        nStairsDone += 1
        print 'finished staircase ' + thisStair.extraInfo['label']
        writeStair(thisStair, filePath)
        if nStairsDone == 1: df = dfStair(thisStair, projName, expNum, targEyeR)
        else: df = pd.concat([df, dfStair(thisStair, projName, expNum, targEyeR)])
        # Recording the data to a csv file:
        df.to_csv(outFileName, index=False, columns=dataCols)
        print 'wrote to ' + outFileName
        completedStairs.append(thisStair)
        print '-------------------------------------------------'

    else:
        nTrialsDone  += 1
        if trialNfb:
            trialNfbText.text = str(nTrialsDone )
        trialNstr = '#' + str(nTrialsDone )

        ## Setting up trial variables

        # mc mask:
        mcSz = thisStair.extraInfo['mcSz']
        mcSf = thisStair.extraInfo['mcSf']
        mcBv = thisStair.extraInfo['mcBv']
        mcBsf = thisStair.extraInfo['mcBsf']
        mcContr = thisStair.extraInfo['mcContr']
        #thisMaskOri = np.rint(np.random.rand(1) * 360)
        #print 'mask ori=' + str(thisMaskOri)
        if mcBv <= 0.01:
            thisMaskFrame = np.random.randint(60)

        # target:
        targSz = thisStair.extraInfo['targSz']
        targSf = thisStair.extraInfo['targSf']
        targYoff = thisStair.extraInfo['targYoff']
        targV = thisStair.extraInfo['targV']

        # random target features: orientation and location:
        targOri1 = thisStair.extraInfo['targOri1']
        targOri2 = thisStair.extraInfo['targOri2']
        allTargOris = np.array([targOri1, targOri2])
        thisTargOri = allTargOris[np.random.randint(2)]
        targXoff1 = thisStair.extraInfo['targXoff1']
        targXoff2 = thisStair.extraInfo['targXoff2']
        allTargXoffs = np.array([targXoff1, targXoff2])
        thisTargXoff = allTargXoffs[np.random.randint(2)]
        targDir = np.random.randint(2) * 2 - 1 # either -1 (up) or 1 (down)

        # setting up the target with the above characteristics:
        targGab.size = targSz
        targGab.sf = targSf
        targGab.ori = thisTargOri + 90 # since, by default, 0 deg results in vert

        # dealing with which eye the targ/mask are presented to:
        if domTest: # else is assigned through GUI
            targEyeR = thisStair.extraInfo['targEyeR']
        maskPos = [winOffX-(2*winOffX*targEyeR), winOffY]
        targGab.pos = [-winOffX+(2*winOffX*targEyeR), winOffY] + \
                      np.array([thisTargXoff, targYoff])

        # temporal variables:
        targTtot = thisStair.extraInfo['targTtot']
        targTpeak = thisStair.extraInfo['targTpeak']
        targTstart = targTpeak-(targTtot/2)
        targTend = targTpeak+(targTtot/2)
        trialT = thisStair.extraInfo['trialT'] # -win.monitorFramePeriod*0.75
        
        print 'TRIAL' + '\t' + 'CONTRAST' + '\t\t' + 'mcBv' + '\t' + 'targT' + \
              '\t' + 'targOri' + '\t' + 'targDir' + '\t' + 'targXoff'
        print trialNstr + '\t' + contrStr + '\t' + str(mcBv) + '\t' + str(targTpeak) + \
              '\t' + str(thisTargOri) + '\t' + str(targDir) + '\t' + str(thisTargXoff)

        # view setup: fade, gap, and fixation cross
        fixCross = thisStair.extraInfo['fixCross']
        periFade = thisStair.extraInfo['mcPeriFade']
        periGap = np.int(mcSz / 2 - periFade)
        #print 'periFade=' + str(periFade) + '; periGap=' + str(periGap)

        nFrames = 60 # number of frames per sequence
        
        # creating an empty matrix for keeping the behavioural responses:
        behRespTrial = []
            
        # initiating the mc gratings:
        mcV = 0
        if mcBv > 0.01:
            grt = np.load(precompiledDir + os.sep + 'mc_' + str(mcV) +
                    '_sf' + str(mcSf) + '_bsf' + str(mcBsf) + '_bv' + str(mcBv) + 
                    '_sz' + str(mcSz) + '.npy')
        else:
            grt = np.load(precompiledDir + os.sep + 'mc_' + str(mcV) +
                    '_sf' + str(mcSf) + '_bsf' + str(mcBsf) + '_bv' + str(9.6) + 
                    '_sz' + str(mcSz) + '.npy')

        # creating a mask, which is fixed for a given trial:
        mcPeriMask = periMask(periGap, periFade)

        #------prepare to start routine "trial"-------
        t = 0
        trialClock.reset()  # clock 
        frameN = -1

        # anchors:
        elStopped = False
        keyPause = False
        targRespGiven = False
        behRespRecorded = False

        # update component parameters for each repeat
        key_arrow = event.BuilderKeyResponse()  # create an object of type keyresponse
        key_arrow.status = NOT_STARTED
        # keep track of which components have finished
        trialComponents = []
        trialComponents.append(winL)
        trialComponents.append(winR)
        trialComponents.append(pauseTextL)
        trialComponents.append(pauseTextR)
        trialComponents.append(fixL)
        trialComponents.append(fixR)
        trialComponents.append(key_arrow)
        trialComponents.append(ISI)
        for thisComponent in trialComponents:
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        
        #-------Start Routine "trial"-------
        continueRoutine = True
        while continueRoutine:
            # get current time
            t = trialClock.getTime()
            frameN = frameN + 1 # number of completed frames (0 is the first frame)
            # update/draw components on each frame
            
            # *winL* updates
            if winL.status == NOT_STARTED:
                # keep track of start time/frame for later
                winL.tStart = t  # underestimates by a little under one frame
                winL.frameNStart = frameN  # exact frame index
                winL.setAutoDraw(True)
                winL.status = STARTED
            
            # *winR* updates
            if winR.status == NOT_STARTED:
                # keep track of start time/frame for later
                winR.tStart = t  # underestimates by a little under one frame
                winR.frameNStart = frameN  # exact frame index
                winR.setAutoDraw(True)
                winR.status = STARTED

            if trialNfb:
                trialNfbText.draw()

            # mcMask and targ presentation:
            if t < trialT:
                # mcMask:
                if mcBv > 0.01:
                    mcMask = visual.GratingStim(win, tex=grt[:,:,frameN%nFrames], 
                        size=(grtSize,grtSize), pos=maskPos, interpolate=False, mask=mcPeriMask)
                else:
                    mcMask = visual.GratingStim(win, tex=grt[:,:,thisMaskFrame], 
                        size=(grtSize,grtSize), pos=maskPos, interpolate=False, mask=mcPeriMask)
                    #ori=thisMaskOri)
                mcMask.contrast = mcContr
                mcMask.draw()
                # drawing the fixation cross, if any:
                if fixCross:
                    fixL.draw()
                    fixR.draw()
                # target presentation:
                if t > targTstart and t < targTpeak:
                    targGab.opacity = sigmoidMod((t-targTstart)*2/targTtot)*10**thisContr
                elif t > targTpeak and t < targTend:
                    targGab.opacity = sigmoidMod((targTend-t)*2/targTtot)*10**thisContr
                else:
                    targGab.opacity = 0
                if t > targTstart and t < targTend:
                    targGab.phase = targGab.phase + (targDir * targV / 60)
                targGab.draw()

            # *key_arrow* updates for target reponses:
            if key_arrow.status == NOT_STARTED and t >= trialT: # record key strokes only after trial end
                # keep track of start time/frame for later
                key_arrow.tStart = t  # underestimates by a little under one frame
                key_arrow.frameNStart = frameN  # exact frame index
                key_arrow.status = STARTED
                # keyboard checking is just starting
                key_arrow.clock.reset()  # now t=0
                event.clearEvents(eventType='keyboard')
                kb_device.clearEvents()
            # registering response at the end of the trial
            if key_arrow.status == STARTED and t > trialT:
                theseKeys = event.getKeys(keyList=['left','right','down'])
                if len(theseKeys) > 0:
                    if 'left' in theseKeys:
                        print 'response: left'
                        behRespTrial = targXoff1 # the first number is negative
                        fdbTextL.text = fdbStrL
                        fdbTextR.text = fdbStrL
                        targRespGiven = True
                    elif 'right' in theseKeys:
                        print 'response: right'
                        behRespTrial = targXoff2 # the second number is positive
                        fdbTextL.text = fdbStrR
                        fdbTextR.text = fdbStrR
                        targRespGiven = True
                    elif 'down' in theseKeys:
                        print 'response: down'
                        behRespTrial = 0
                        fdbTextL.text = fdbStr0
                        fdbTextR.text = fdbStr0
                        targRespGiven = True
                    if targRespGiven: # this is overwritten every time any key is pressed
                        rt = t - targTstart
                        if behRespTrial == thisTargXoff: corrResp = 1 # correct dir resp
                        else: corrResp = 0 # this is the result from both 0 and wrong dir resp
                        if thisContr <= -2: corrResp = 0

            if t > trialT and not keyPause:
                fdbTextL.draw()
                fdbTextR.draw()

            # pause text and data exporting
            if targRespGiven and not keyPause and t>trialT:
                pauseTextL.draw()
                pauseTextR.draw()
                if 'space' in event.getKeys(keyList=['space']):
                    keyPause = True
                    # only update the response upon pressing 'space':
                    fdbTextL.text = fdbStrQn
                    fdbTextR.text = fdbStrQn
                    thisStair.addResponse(corrResp)
                    thisStair.addOtherData('rt',rt)
                    thisStair.addOtherData('targOri',thisTargOri)
                    thisStair.addOtherData('targXoff',thisTargXoff)
                    stairs.append(thisStair)
                    #print 'spacebar pressed'
                    pauseTextL.setAutoDraw(False)
                    pauseTextR.setAutoDraw(False)

            # *ISI* period
            if ISI.status == NOT_STARTED and t>=trialT and keyPause:
                # keep track of start time/frame for later
                ISI.tStart = t  # underestimates by a little under one frame
                ISI.frameNStart = frameN  # exact frame index
                fixL.setAutoDraw(True)
                fixR.setAutoDraw(True)
                ISI.start(ISIduration)
            #one frame should pass before updating params and completing
            elif ISI.status == STARTED and t >= (ISI.tStart + ISIduration): 
                fixL.setAutoDraw(False)
                fixR.setAutoDraw(False)
                ISI.complete() #finish the static period
                continueRoutine = False
            
            # check if all components have finished
            # a component has requested a forced-end of Routine:
            if not continueRoutine: 
                # if we abort early the non-slip timer needs reset:
                routineTimer.reset() 
                break
            # will revert to True if at least one component still running
            continueRoutine = False  
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and \
                        thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # check for quit (the Esc key)
            if endExpNow or event.getKeys(keyList=["escape"]):
                print np.shape(behResp)
                core.quit()
            
            # refresh the screen
            # don't flip if this routine is over or we'll get a blank screen
            if continueRoutine:  
                win.flip()
            else: # this Routine was not non-slip safe so reset non-slip timer
                routineTimer.reset()
        
        #-------Ending Routine "trial"-------
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)

    # thisExp.nextEntry()

# summary of completed stairs:
stairList, eyeList, dfThresh = [], [], []
for thisStair in completedStairs:
    stairList.append(thisStair.extraInfo['label'])
    dfThresh.append(thisStair.reversalIntensities)
    print 'meanRev6 = %.3f' %(np.average(thisStair.reversalIntensities[-6:]))
    #print thisStair.reversalIntensities
dfThresh = pd.DataFrame.from_records(dfThresh, index=stairList)
print dfThresh
dfThresh.to_csv(threshFileName, index=stairList, header=False)

if expNum == 'dom': # at the end of dom test, display averaged thresh's for each eye:
    print dfThreshDomRaw(dfThresh)
    
print "finished the experiment"

win.close()
core.quit()
