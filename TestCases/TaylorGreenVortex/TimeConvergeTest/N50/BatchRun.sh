#! /bin/bash

cd DT0.005
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.0025
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.00125
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.000625
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.0003125
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.00015625
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.000078125
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

