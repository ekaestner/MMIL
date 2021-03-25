MATLAB Compiler

1. Prerequisites for Deployment 

. Verify the MATLAB Compiler Runtime (MCR) is installed and ensure you    
  have installed version 8.2 (R2013b).   

. If the MCR is not installed, do the following:
  (1) enter
  
      >>mcrinstaller
      
      at MATLAB prompt. The MCRINSTALLER command displays the 
      location of the MCR Installer.

  (2) run the MCR Installer.

Or download the Linux 64-bit version of the MCR for R2013b 
from the MathWorks Web site by navigating to

   http://www.mathworks.com/products/compiler/mcr/index.html
   
   
For more information about the MCR and the MCR Installer, see 
Distribution to End Users in the MATLAB Compiler documentation  
in the MathWorks Documentation Center.    


2. Files to Deploy and Package

Files to package for Standalone 
================================
-muxrecon 
-run_muxrecon.sh (shell script for temporarily setting environment variables and 
                  executing the application)
   -to run the shell script, type
   
       ./run_muxrecon.sh <mcr_directory> <argument_list>
       
    at Linux or Mac command prompt. <mcr_directory> is the directory 
    where version 8.2 of MCR is installed or the directory where 
    MATLAB is installed on the machine. <argument_list> is all the 
    arguments you want to pass to your application. For example, 

    If you have version 8.2 of the MCR installed in 
    /mathworks/home/application/v82, run the shell script as:
    
       ./run_muxrecon.sh /mathworks/home/application/v82
       
    If you have MATLAB installed in /mathworks/devel/application/matlab, 
    run the shell script as:
    
       ./run_muxrecon.sh /mathworks/devel/application/matlab
-MCRInstaller.zip
   -if end users are unable to download the MCR using the above  
    link, include it when building your component by clicking 
    the "Add MCR" link in the Deployment Tool
-This readme file 

3. Definitions

For information on deployment terminology, go to 
http://www.mathworks.com/help. Select MATLAB Compiler >   
Getting Started > About Application Deployment > 
Application Deployment Terms in the MathWorks Documentation 
Center.


4. Appendix 

A. Linux x86-64 systems:
   On the target machine, add the MCR directory to the environment variable 
   LD_LIBRARY_PATH by issuing the following commands:

        NOTE: <mcr_root> is the directory where MCR is installed
              on the target machine.         

            setenv LD_LIBRARY_PATH
                $LD_LIBRARY_PATH:
                <mcr_root>/v82/runtime/glnxa64:
                <mcr_root>/v82/bin/glnxa64:
                <mcr_root>/v82/sys/os/glnxa64
            setenv XAPPLRESDIR <mcr_root>/v82/X11/app-defaults

   For more detail information about setting MCR paths, see Distribution to End Users in 
   the MATLAB Compiler documentation in the MathWorks Documentation Center.


     
        NOTE: To make these changes persistent after logout on Linux 
              or Mac machines, modify the .cshrc file to include this  
              setenv command.
        NOTE: The environment variable syntax utilizes forward 
              slashes (/), delimited by colons (:).  
        NOTE: When deploying standalone applications, it is possible 
              to run the shell script file run_muxrecon.sh 
              instead of setting environment variables. See 
              section 2 "Files to Deploy and Package".    






