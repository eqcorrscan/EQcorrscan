@ECHO OFF
call "%WIN_SDK_ROOT%\%WINDOWS_SDK_VERSION%\Setup\WindowsSdkVer.exe" -q -version:%WINDOWS_SDK_VERSION%
call "%WIN_SDK_ROOT%\%WINDOWS_SDK_VERSION%\Bin\SetEnv.cmd"

ECHO Downloading DLL files
IF %PYTHON_ARCH% == 64 (
    call curl.exe -sS -o fftw-3.3.5-dll64.zip "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.5-dll64.zip"
    SET MACHINE=X64
    SET FFTW_DLL_FILENAME=fftw-3.3.5-dll64.zip
) ELSE (
    call curl.exe -sS -o fftw-3.3.5-dll32.zip "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.5-dll32.zip"
    SET MACHINE=X86
    SET FFTW_DLL_FILENAME=fftw-3.3.5-dll32.zip
)
ECHO Extracting DLLs from %FFTW_DLL_FILENAME%
call 7z.exe e %FFTW_DLL_FILENAME% -oeqcorrscan\utils\lib *.dll
call 7z.exe e %FFTW_DLL_FILENAME% -oeqcorrscan\utils\lib *.def
ECHO Generating def files
call lib /machine:%MACHINE% /def:eqcorrscan\utils\lib\libfftw3-3.def
call lib /machine:%MACHINE% /def:eqcorrscan\utils\lib\libfftw3f-3.def
call lib /machine:%MACHINE% /def:eqcorrscan\utils\lib\libfftw3l-3.def