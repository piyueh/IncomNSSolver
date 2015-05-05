#! /bin/bash

cd "Tol1E-15"
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m ../meshSetting.dat -d ../initData.dat -c compParams.dat &> log.txt &
cd ..

cd "Tol1E-13"
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m ../meshSetting.dat -d ../initData.dat -c compParams.dat &> log.txt &
cd ..

cd "Tol1E-11"
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m ../meshSetting.dat -d ../initData.dat -c compParams.dat &> log.txt &
cd ..

cd "Tol1E-9"
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m ../meshSetting.dat -d ../initData.dat -c compParams.dat &> log.txt &
cd ..

cd "Tol1E-7"
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m ../meshSetting.dat -d ../initData.dat -c compParams.dat &> log.txt &
cd ..

cd "Tol1E-5"
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m ../meshSetting.dat -d ../initData.dat -c compParams.dat &> log.txt &
cd ..

cd "Tol1E-3"
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m ../meshSetting.dat -d ../initData.dat -c compParams.dat &> log.txt &
cd ..

