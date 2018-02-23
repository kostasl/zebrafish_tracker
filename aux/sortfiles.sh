##Batch script to move  files in correct directory structure
mkdir DryFed
mkdir LiveFed
mkdir NotFed
mkdir DryFed/Empty/
mkdir DryFed/Live/
mkdir NotFed/Live/
mkdir NotFed/Empty/
mkdir LiveFed/Empty/
mkdir LiveFed/Live/

cp *DryFed*Empty*.* DryFed/Empty/ -v
cp *DryFed*Live*.* DryFed/Live/ -v
cp *DryFed*Roti*.* DryFed/Live/ -v
cp *LiveFed*Empty*.* LiveFed/Empty/ -v
cp *LiveFed*Roti*.* LiveFed/Live/ -v
cp *LiveFed*Live*.* LiveFed/Live/ -v
cp *NotFed*Roti*.* NotFed/Live/ -v
cp *NotFed*Live*.* NotFed/Live/ -v
cp *NotFed*Empty*.* NotFed/Empty/ -v
