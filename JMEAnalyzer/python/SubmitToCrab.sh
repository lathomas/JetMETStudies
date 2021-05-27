#!/bin/bash
thedataset=$1
campaigntag=`echo $thedataset | sed -e 's/.*\/\(.*\)\/MINIAOD.*/\1/'`
campaigntag=`echo $campaigntag | sed -e 's/[-]/\_/g'`
thedataset=`echo $1 | sed -e 's/[-]/HYPHEN/g'`
thedataset=`echo $thedataset | sed -e 's/\//SLASH/g'`
thedataset=`echo $thedataset | sed -e 's/\_/UNDERSCORE/g'`
njobs=$3
skim=$4


therequestname=$2


ISMC="True"
echo "$thedataset" | grep "AODSIM" && ISMC="True" || ISMC="False"
ISFS="False"
echo "$thedataset" | grep "Fast" && ISFS="True" || ISFS="False"



RUNERA=$5
echo "runera is $RUNERA"

cp crab_template.py crab_template_temp.py

sed -ie "s/THEDATASET/$thedataset/g" crab_template_temp.py

sed -ie 's/HYPHEN/-/g' crab_template_temp.py
sed -ie 's/SLASH/\//g' crab_template_temp.py
sed -ie 's/UNDERSCORE/\_/g' crab_template_temp.py

sed -ie "s/NJOBS/$njobs/g" crab_template_temp.py
sed -ie "s/THESKIM/$skim/g" crab_template_temp.py
sed -ie "s/RUNERA/$RUNERA/g" crab_template_temp.py
sed -ie "s/THEREQUESTNAME/$therequestname/g" crab_template_temp.py
sed -ie "s/CAMPAIGN/$campaigntag/g" crab_template_temp.py



cp JMEanalysis_forcrab.py theconfig.py 
sed -ie "s/THESKIM/$skim/g" theconfig.py
sed -ie "s/MCBOOL/$ISMC/g" theconfig.py
sed -ie "s/FSBOOL/$ISFS/g" theconfig.py
sed -ie "s/THERUNERA/$RUNERA/g" theconfig.py

#cp theconfig.py  theconfigbu.py

crab submit crab_template_temp.py 
rm theconfig.py crab_template_temp.py



