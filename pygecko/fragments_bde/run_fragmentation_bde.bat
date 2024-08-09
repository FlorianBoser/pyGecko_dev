@echo off

:: Activate the BDE environment
call conda activate pyGecko_X_BDE

:: Set the PYTHONPATH to include the directory where the BDE fragments module is located
set PYTHONPATH=%PYTHONPATH%;C:\Users\flori\anaconda3\envs\pyGecko\Lib\site-packages\pygecko\fragments_bde

:: Debug PYTHONPATH
echo PYTHONPATH is: %PYTHONPATH%

:: Run the BDE service script with the passed arguments
python "C:\Users\flori\anaconda3\envs\pyGecko\Lib\site-packages\pygecko\fragments_bde\fragments_bde.py" %*