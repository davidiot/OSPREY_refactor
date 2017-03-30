/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.TreeMap;

//This class stores parameters from input files
//Control package classes will probably each want a param set
//which will control its choice of algorithms, and will be used to initialize other classes 
//like EnergyMatrix and A* according to user settings


import edu.duke.cs.osprey.handlempi.MPIMaster;
import edu.duke.cs.osprey.tools.StringParsing;

/**
 * Handles reading in and managing parameter/value pairs from the input configuration files
 */
public class ParamSet implements Serializable {
	
	private static final long serialVersionUID = 4364963601242324780L;

	private TreeMap<String,String> params = new TreeMap<>();//map parameter/value pairs
        //parameter names will be stored as all upper-case, to avoid confusion
	
        private static final String defaultParamFile = "defaults.cfg";
        private TreeMap<String,String> defaultParams = new TreeMap<>();
        
        private TreeMap<String,String> wildcardDefaults = new TreeMap<>();
        //Handles default params of the form "BLABLA* 1" (would map BLABLA -> 1)
        //then for example if we needed a default for parameter BLABLABLA, it would return 1
        
        private boolean isVerbose;
        
	//constructor
	ParamSet(){
       isVerbose = true;     
	}
	
	public void setVerbosity(boolean val) {
		isVerbose = val;
	}
	
	public boolean isVerbose() {
		return isVerbose;
	}
	
	//Reads in all parameter pairs from the file fName and updates the params
	public void addParamsFromFile(String fName){
            loadParams(fName, params);
        }
        
        public void addDefaultParams(){
            String defaultFilePath = EnvironmentVars.getDataDir() + defaultParamFile;
            loadParams(defaultFilePath, defaultParams);
            
            for(String param : defaultParams.keySet()){
                if(param.endsWith("*")){
                    String wildcard = param.substring(0, param.length()-1);
                    wildcardDefaults.put( wildcard, defaultParams.get(param) );
                }
            }
        }
	
        
	private static void loadParams(String fName, TreeMap<String,String> paramMap) {
		// load all parameters for cfg file fName into the map paramMap
		
		// First attempt to open and read the config file
		try (BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(fName)))) {
			
			// read line-by-line
			String line;
			while ((line = in.readLine()) != null) {
				
				// strip comments
				int commentStartPos = line.indexOf('%');
				if (commentStartPos >= 0) {
					line = line.substring(0, commentStartPos);
				}
				
				// skip blank lines
				if (line.isEmpty()) {
					continue;
				}
				
				// parse the param
				String paramName = StringParsing.getToken(line, 1).trim();
				paramName = paramName.toUpperCase();
				if (paramMap.containsKey(paramName)) {
				
					// duplicate param
					throw new RuntimeException("ERROR: parameter " + paramName + " already read");
					
				} else {
					
					// new param
					String paramVal = line.substring(paramName.length() + 1);
					paramMap.put(paramName, paramVal);
				}
			}
			
		} catch (FileNotFoundException ex) {
			throw new RuntimeException("ERROR: Couldn't find configuration file " + fName);
		} catch (Exception ex) {
			throw new Error("ERROR: An error occurred reading configuration file " + fName, ex);
		}
	}

	
	//Sets the value of parameter paramName to newValue
	public void setValue(String paramName, String newValue){
		paramName = paramName.toUpperCase();
		params.put(paramName, newValue);
	}

        
        
        //Methods to get parameter values
        //We first try to get a value from params (loaded from user's cfg files)
        //If there is none, we try to get a value from defaultParams (loaded from default cfg)
        //If there is none there either, this is an error
        
	public String getValue(String paramName){
            
            paramName = paramName.toUpperCase();
            String val = params.get(paramName);
            
            if(val==null){
                val = defaultParams.get(paramName);
                
                if(val==null){//See if this is a wildcard param
                    for(String wildcard : wildcardDefaults.keySet()){
                        if(paramName.startsWith(wildcard)){
                            val = wildcardDefaults.get(wildcard);
                            break;
                        }
                    }
                }
                
                if(val==null)//still null
                    throw new RuntimeException("ERROR: Parameter "+paramName+" not found");
                
                if(isVerbose)
                	MPIMaster.printIfMaster("Parameter "+paramName+" not set. Using default value "+val);
            }
            else {
            	if(isVerbose)
            		MPIMaster.printIfMaster("Parameter "+paramName+" set to "+val);
            }
            
            return val.trim();
	}
        
        
        //The following methods return a parameter expected to be a certain non-string type
        //It's an error if they're not that type
        public int getInt(String paramName){
            String val = getValue(paramName);
            try {
                return Integer.valueOf(val);
            }
            catch(NumberFormatException e){
                throw new RuntimeException("ERROR: Value "+val+" for parameter "+paramName+" can't be parsed as an integer");
            }
        }
        
        public boolean getBool(String paramName){
            String val = getValue(paramName);
            try {
                return Boolean.valueOf(val);
            }
            catch(NumberFormatException e){
                throw new RuntimeException("ERROR: Value "+val+" for parameter "+paramName+" can't be parsed as a boolean");
            }
        }
        
        public double getDouble(String paramName){
            String val = getValue(paramName);
            try {
                return Double.valueOf(val);
            }
            catch(NumberFormatException e){
                throw new RuntimeException("ERROR: Value "+val+" for parameter "+paramName+" can't be parsed as a double");
            }
        }
        
        
        public String getRunSpecificFileName(String paramName, String suffix){
            //This is for a special kind of parameter whose value is 
            //the name of a run-specific file.  So by default the value is runName + suffix
            paramName = paramName.toUpperCase();
            
            if(params.containsKey(paramName))
                return params.get(paramName);
            else
                return getValue("RUNNAME") + suffix;
        }
        
        
        
        
        //searching for a parameter
        public ArrayList<String> searchParams(String searchTerm){
            //return any parameter names that include searchTerm
            searchTerm = searchTerm.toUpperCase();//all parameters are upper-case
            
            ArrayList<String> matches = new ArrayList<>();//parameters matching the searchTerm
            
            for(String paramName : params.keySet()){
                if(paramName.contains(searchTerm)){
                    matches.add(paramName);
                }
            }
            
            return matches;
        }
        
        

      //Methods to get parameter values
        //Default methods can be used if not set; if null default then return an error
        
	public String getValue(String paramName, String defaultVal){
            
            paramName = paramName.toUpperCase();
            String val = params.get(paramName);
            
            if(val==null){
                if(defaultVal==null)//no default...must be set
                    throw new RuntimeException("ERROR: Parameter "+paramName+" not found");
                
                val = defaultVal;
                if(isVerbose) 
                	MPIMaster.printIfMaster("Parameter "+paramName+" not set. Using default value "+defaultVal);
            }
            else {
            	if(isVerbose)
            		MPIMaster.printIfMaster("Parameter "+paramName+" set to "+val);
            }
            
            return val.trim();
	}
        
        
        //The following methods return a parameter expected to be a certain non-string type
        //It's an error if they're not that type
        public int getInt(String paramName, int defaultVal){
            String val = getValue(paramName,String.valueOf(defaultVal));
            try {
                return Integer.valueOf(val);
            }
            catch(NumberFormatException e){
                throw new RuntimeException("ERROR: Value "+val+" for parameter "+paramName+" can't be parsed as an integer");
            }
        }
        
        public boolean getBool(String paramName, boolean defaultVal){
            String val = getValue(paramName,String.valueOf(defaultVal));
            try {
                return Boolean.valueOf(val);
            }
            catch(NumberFormatException e){
                throw new RuntimeException("ERROR: Value "+val+" for parameter "+paramName+" can't be parsed as a boolean");
            }
        }
        
        public double getDouble(String paramName, double defaultVal){
            String val = getValue(paramName,String.valueOf(defaultVal));
            try {
                return Double.valueOf(val);
            }
            catch(NumberFormatException e){
                throw new RuntimeException("ERROR: Value "+val+" for parameter "+paramName+" can't be parsed as a double");
            }
        }
        
}
