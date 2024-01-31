#!/usr/bin/perl -w
# Copy renamed anainc 30 day statistics plots to montly archive directory
# Version 2.0 / MKa 14.4.2022
# - v 2.0 : cycle subdivison added
# - 20221104 : version for MNWC prod (paths)


if ($ARGV[2]){
    $YYYY = $ARGV[0];
    $MM = $ARGV[1];
    $HH = $ARGV[2];
} else {
    print "\nPlease give three arguments : YYYY  MM  HH \nAborting\n\n";
    exit;
};

$tStamp = $YYYY.$MM;
$pathFrom = "/data/hirlam2/Python_MNWC/".$HH."/aninc/";
$pathTo = "/fmi/data/fminwp.fmi.fi/Harmonie-pyth-MNWC/aninc_monthly/".$YYYY."/".$MM."/".$HH."/";;

print "$pathFrom \n";
print "$tStamp \n";
print "$pathTo \n";

chdir($pathFrom);
@listaus = `ls -1 *_average*`;

@args = ("mkdir","-p","$pathTo");       #    * create YEAR/MONTH/HOUR directory ...
if (!-e $pathTo) {system(@args)};       #      ... if it does not exist

foreach $fn (@listaus){
    chomp ($fn);
    $fnsrc = $fn;
    $fnsrc =~ /(\w*).png/;
    $fnew = $1."_".$tStamp.".png";   # fname.png => fname_tStamp.png
    print "$fnsrc  $1  $fnew \n";
    `cp -p $fn $pathTo\/$fnew`;
};

exit;


