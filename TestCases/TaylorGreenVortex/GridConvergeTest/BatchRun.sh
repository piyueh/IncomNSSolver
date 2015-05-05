#! /bin/bash

cd N25
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m meshSetting.dat -d initData.dat -c ../compParams.dat &> log.txt &
cd ..

cd N51
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m meshSetting.dat -d initData.dat -c ../compParams.dat &> log.txt &
cd ..

cd N101
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m meshSetting.dat -d initData.dat -c ../compParams.dat &> log.txt &
cd ..

cd N151
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m meshSetting.dat -d initData.dat -c ../compParams.dat &> log.txt &
cd ..

cd N201
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m meshSetting.dat -d initData.dat -c ../compParams.dat &> log.txt &
cd ..

cd N251
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m meshSetting.dat -d initData.dat -c ../compParams.dat &> log.txt &
cd ..

cd N301
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m meshSetting.dat -d initData.dat -c ../compParams.dat &> log.txt &
cd ..

cd N351
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m meshSetting.dat -d initData.dat -c ../compParams.dat &> log.txt &
cd ..

cd N401
rm *.txt
nohup IncomNSSolver -f ../fluid.dat -m meshSetting.dat -d initData.dat -c ../compParams.dat &> log.txt &
cd ..

