#! /bin/bash

cd Re100
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

cd Re1000
rm *.txt
nohup IncomNSSolver -f fluid.dat -m meshSetting.dat -d initData.dat -c compParams.dat &> log.txt &
cd ..

