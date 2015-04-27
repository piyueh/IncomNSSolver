#! /bin/bash

cd DT0.0008
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d ../N401Data.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.0004
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d ../N401Data.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.0002
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d ../N401Data.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.0001
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d ../N401Data.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.00005
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d ../N401Data.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.000025
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d ../N401Data.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.0000125
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d ../N401Data.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.00000625
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d ../N401Data.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.000003125
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d ../N401Data.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.0000015625
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d ../N401Data.dat -c compParams.dat &> log.txt &
cd ..

cd DT0.00000078125
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d ../N401Data.dat -c compParams.dat &> log.txt &
cd ..


