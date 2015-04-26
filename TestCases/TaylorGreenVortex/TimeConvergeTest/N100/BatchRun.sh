#! /bin/bash

cd DT0.01
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.005
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.0025
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.00125
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.000625
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.0003125
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.00015625
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.000078125
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

