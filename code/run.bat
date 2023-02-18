Setlocal EnableDelayedExpansion
echo off
g++ main.cpp -o main -fopenmp
for /r %%i in (/test/*.ppm) do (
   set boba=%%i
   set in=!boba:\code\=\code\test\!
   set out=!in:.ppm=_out.ppm!
   set out=!out:\test\=\result\!
   echo !in!
   main.exe 8 !in! !out! 0.0
)
PAUSE