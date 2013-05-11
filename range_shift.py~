import sys, time, string, os, math, re
import ConfigParser
import traceback

import numpy

from sys import path
from string import split
from math import sqrt
from numpy import *

import copy
import shutil
import subprocess

from current_vectors import *

version = '2013_0408'

options = {}
try:
    ####### USER SETTINGS ##############################
    circuitscapeDir = "c:\\temp" # Directory where Circuitscape directory is placed
    arrowOptions[writeTotalCurrent] = False # Write total current leaving each pixel
    arrowOptions[writeResultant] = True # Write vector magnitudes (will be less than total current)
    arrowOptions[writeAllDirections] = False # Write current leaving pixels from N, NE, E, etc..
    arrowOptions[writeEdgeMaps] = False # in Progress
    arrowOptions[deleteTempFiles] = True # Delete everything except standard current map (saved for debug)
    arrowOptions[writeArcVectors] = False
    raiseResistPower = 1
    
    if len(sys.argv) < 2: #Manual inputs
        resistInput = 'c:\\dropbox\\working\\arrows\\nocost.asc' # raw resistance map (no zeros allowed, NoData = -9999)
        raiseResistPower = 1
        # resistInput = 'c:\\arrows\\r3.asc' # raw resistance map (no zeros allowed, NoData = -9999)
        rangeDir = 'C:\\dropbox\\working\\a1332' # All range maps in this directory, NoData = -9999
        # rangeDir = 'c:\\arrows\\testranges' # All range maps in this directory, NoData = -9999
        scratchDir = 'c:\\dropbox\\working\\arrows\\scratch'
        outDir = 'c:\\temp2'
        ####################################################
        fileList = os.listdir(rangeDir)
        settings = {}
        settings['resistInput'] = resistInput  
        settings['rangeDir'] = rangeDir
        settings['scratchDir'] = scratchDir
        settings['outDir'] = outDir

    else:
        resistInput = str(sys.argv[1])
        # raw resistance map (no zeros allowed, NoData = -9999)
        rangeDir = str(sys.argv[2]) # All range maps in this directory, NoData = -9999
        scratchDir = str(sys.argv[3])
        rootOut = str(sys.argv[4])
        taxa = str(sys.argv[5])
        gcm = str(sys.argv[6])
        resistance = str(sys.argv[7])
        settings = str(sys.argv)
        outDir = rootOut + "/" + resistance + "/" + taxa + "/" + gcm  + "/a2/2071-2100"

        fileListRaw = os.listdir(rangeDir)
        fileList = list()
        for filename in fileListRaw:
            match = re.match(r'[ABM][1-9]*.asc$', filename)
            if match is not None:
                fileList.append(filename)
    print 'Creating directories'
    if not os.path.exists(scratchDir):
        os.mkdir(scratchDir)
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    logDir = os.path.join(scratchDir,'log')    
    if not os.path.exists(logDir):
        os.mkdir(logDir)
   
    logFilePath = create_log_file(logDir, settings)    
    lprint(logFilePath, 'Version = ' + version + '\n')
    
    for fileName in fileList:   
        baseFN, ext = os.path.splitext(fileName)
        if ext != '.asc':
            continue
        lprint(logFilePath, 'Processing ' + fileName)
        state = {}
        baseResistFN, ext = os.path.splitext(resistInput)
        projectionFile = baseResistFN + '.prj'
        rangeRaster = os.path.join(rangeDir,fileName)
        state,rangeMap = read_map(state,rangeRaster)
        contractResistFile = os.path.join(scratchDir, baseFN + '_contractResist.asc')
        expandResistFile = os.path.join(scratchDir, baseFN + '_expandResist.asc')
        sourceFile = os.path.join(scratchDir, baseFN + '_sources.asc')
        targetFile = os.path.join(scratchDir, baseFN + '_targets.asc')
        expandGroundFile = os.path.join(scratchDir, baseFN + '_expand_grounds.asc')
        contractGroundFile = os.path.join(scratchDir, baseFN + '_contract_grounds.asc')
        expandConfigFile = os.path.join(scratchDir, baseFN + '_expand.ini')
        contractConfigFile = os.path.join(scratchDir, baseFN + '_contract.ini')
        state, resisMap = read_map(state, resistInput)
        
        if raiseResistPower !=1:
            lprint(logFilePath,'*************************************')
            lprint(logFilePath,'Raising resistance by power of '+ str(raiseResistPower))
            resisMap = where(resisMap > 0, power(resisMap,raiseResistPower), resisMap)
            
            # File naming below is just for naming output files in current vectors script
            baseResistFN, ext = os.path.splitext(resistInput)
            resistPowerFile = baseResistFN + 'e' + str(raiseResistPower) +  ext 
            resistInput = resistPowerFile      
                
        contractResist = where(rangeMap==1,-9999,resisMap)
        contractResist = where(rangeMap==3,-9999,contractResist)
        contractResist = where(rangeMap==-9999,-9999,contractResist)
        expandResist = where(rangeMap==1,-9999,resisMap)
        expandResist = where(rangeMap==-9999,-9999,expandResist)
        del resisMap
        
        writer(contractResistFile, contractResist, state, False, projectionFile)
        del contractResist
        writer(expandResistFile, expandResist, state, False, projectionFile)
        del expandResist
        sources = where(rangeMap==2,1,0)
        writer(sourceFile, sources, state, False, projectionFile)
        del sources
        targets = where(rangeMap==3,-1,0)
        writer(targetFile, targets, state, False, projectionFile)
        del targets
        contractGrounds = where(rangeMap==4,0,-9999)
        writer(contractGroundFile, contractGrounds, state, False, projectionFile)
        expandGrounds = where(rangeMap==2,0, contractGrounds)
        writer(expandGroundFile, expandGrounds, state, False, projectionFile)
        

        # numExpandCells = len(where(rangeMap==2))
        # print numExpandCells
        del contractGrounds, expandGrounds

        cs_options=set_default_cs_options()
        
        # Contraction
        cs_options['ground_file']= contractGroundFile
        cs_options['habitat_file'] = contractResistFile
        cs_options['source_file']= sourceFile
        cs_options['output_file'] = os.path.join(scratchDir,baseFN+'_contract.out')       
        writeConfigFile(contractConfigFile,cs_options)
        
        # Expansion
        cs_options['ground_file']= expandGroundFile
        cs_options['habitat_file'] = expandResistFile
        cs_options['source_file']= targetFile
        cs_options['output_file'] = os.path.join(scratchDir,baseFN+'_expand.out')
        writeConfigFile(expandConfigFile,cs_options)
        
        # ...call circuitscape twice, once for expanding once for contracting
        CSPATH = get_cs_path(circuitscapeDir)  
        startTime=time.clock()
        lprint(logFilePath,'\nExpansion analysis')
        memFlag1, failFlag1, solverFailFlag1 = call_circuitscape(CSPATH, expandConfigFile, logFilePath)
        if failFlag1 == True:
            lprint(logFilePath, '******************************************')
            lprint(logFilePath, 'Circuitscape failed for expansion analysis')
            lprint(logFilePath, '******************************************\n\n\n')        
        if solverFailFlag1 == True:
            lprint (logFilePath, '********************************************')
            lprint (logFilePath, 'Circuitscape SOLVER failed for expansion analysis. Yarg.')  
            lprint (logFilePath, '******************************************\n\n\n')  
            lprint (logFilePath, 'Settings:')            
            lprint (logFilePath, str(settings))            
            
        lprint(logFilePath,'Contraction analysis')        
        memFlag2, failFlag2, solverFailFlag2 = call_circuitscape(CSPATH, contractConfigFile, logFilePath)
        if failFlag2 == True:
            lprint (logFilePath, '********************************************')
            lprint (logFilePath, 'Circuitscape failed for contraction analysis')  
            lprint (logFilePath, '******************************************\n\n\n')        

        if solverFailFlag2 == True:
            lprint (logFilePath, '********************************************')
            lprint (logFilePath, 'Circuitscape SOLVER failed for contraction analysis. Yarg.')  
            lprint (logFilePath, '******************************************\n\n\n') 
            lprint (logFilePath, 'Settings:')            
            lprint (logFilePath, str(settings))            
        lprint (logFilePath, 'Done with circuitscape calls.')
        startTime2 = elapsed_time(logFilePath,startTime)
        
        # Map vectors
        lprint(logFilePath, 'Calculating vector directions and magnitudes')
        if failFlag1 == False and solverFailFlag1 == False:
            cs_options = readConfigFile(expandConfigFile)
            csOutDir, baseOutputFN = os.path.split(cs_options['output_file'])
            baseFile, ext = os.path.splitext(baseOutputFN)
            expandVoltMapFile = os.path.join(csOutDir, baseFile + '_voltmap.asc')        
            map_current_vectors(expandConfigFile, expandVoltMapFile, arrowOptions, rangeRaster, resistInput, projectionFile, outDir, deleteTempFiles, logFilePath) 
        if failFlag2 == False and solverFailFlag2 == False:
            cs_options = readConfigFile(contractConfigFile)
            csOutDir, baseOutputFN = os.path.split(cs_options['output_file'])
            baseFile, ext = os.path.splitext(baseOutputFN)
            contractVoltMapFile = os.path.join(csOutDir, baseFile + '_voltmap.asc')        
            map_current_vectors(contractConfigFile, contractVoltMapFile, arrowOptions, rangeRaster, resistInput, projectionFile, outDir, deleteTempFiles, logFilePath)
        lprint (logFilePath, '\n-------------------------------------------------------')    
        lprint (logFilePath, 'Done with all operations for ' + fileName +'\n')
        startTime = elapsed_time(logFilePath, startTime)
        # Clean up    
        if deleteTempFiles == True:
            try:
                os.remove(contractResistFile) 
                os.remove(expandResistFile)
                os.remove(sourceFile)
                os.remove(targetFile)
                os.remove(contractGroundFile)
                os.remove(expandGroundFile)
                fileList = os.listdir(scratchDir)  
                for fileName in fileList:   
                    baseFN, ext = os.path.splitext(fileName)
                    if ext == '.prj' and 'curmap' not in baseFN: # Saving standard current maps (when written) for debug
                        os.remove(os.path.join(scratchDir,fileName))
                # shutil.rmtree(scratchDir)
            except:
                pass
    lprint(logFilePath,'\nDone with all operations for directory ' + rangeDir)               

except:
    exit_with_python_error(logFilePath,'range_shift.py')
        
        