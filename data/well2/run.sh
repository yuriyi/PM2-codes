#!/bin/bash
################################################################################
## Test of fkm3d program - elastic/acoustic wave field modeling in horizontally
## stratified media using reflectivity method.
################################################################################

## Medium parameters
# medium.txt

# Receiver parameters
xrec0=0;    yrec0=0;    zrec0=0;    # receiver line beginning
xrecN=1000; yrecN=0;    zrecN=0;    # receiver line end
nrec=51                             # number of receivers
typerec=3                           # receiver type (1=vx, 2=vy, 3=vz, 4=p)

# Source parameters
xsrc=0 ; ysrc=0 ; zsrc=400 ;         # source coordinates
Fx=0.0 ; Fy=0.0 ; Fz=1.0 ;
Mxx=0.0; Mxy=0.0; Mxz=0.0;
Myx=0.0; Myy=0.0; Myz=0.0;
Mzx=0.0; Mzy=0.0; Mzz=0.0;

# Time parameters
nf=101      # number of frequency samples
fim=0.2     # imaginary part for frequencies
nt=2001     # number of time samples
dt=0.0008   # time step
t0=0.1      # shot time
f0=20       # central frequency

# Slowness parameters
nu=1001     # number of slowness samples
du=2.0e-6   # slowness step

#################################################################################################################################################################################
## Auxiliary variables
#################################################################################################################################################################################

df=`python <<END
print 1/(($nt-1)*$dt)
END`
ntan=`python <<END
print 2*($nf-1)
END`
dtan=`python <<END
print 1/($ntan*$df)
END`
fsan=`python <<END
print $ntan*$df
END`
omim=`python <<END
import numpy as np
print 2*np.pi*$fim
END`

dxrec=`python <<END
if $nrec!=0:
    print float($xrecN-$xrec0)/($nrec-1)
else:
    print 0
END`
dyrec=`python <<END
if $nrec!=0:
    print float($yrecN-$yrec0)/($nrec-1)
else:
    print 0
END`
dzrec=`python <<END
if $nrec!=0:
    print float($zrecN-$zrec0)/($nrec-1)
else:
    print 0
END`

vp_ar=($vp)
vs_ar=($vs)
rh_ar=($rh)
om_ar=($om)
Qp_ar=($Qp)
Qs_ar=($Qs)
dd_ar=($dd)

fill-linear-coordinates --ofile receivers.txt --n $nrec --x0 $xrec0 --xN $xrecN --y0 $yrec0 --yN $yrecN --z0 $zrec0 --zN $zrecN --kind $typerec

printf "%f\t%f\t%f\t\t\t\\\\\\ x, y, z\n" {$xsrc,$ysrc,$zsrc} > source.txt
printf "%f\t%f\t%f\t\t\t\\\\\\ Fx, Fy, Fz\n" {$Fx,$Fy,$Fz} >> source.txt
printf "%f\t%f\t%f\t\t\t\\\\\\ Mxx, Mxy, Mxz\n" {$Mxx,$Mxy,$Mxz} >> source.txt
printf "%f\t%f\t%f\t\t\t\\\\\\ Myx, Myy, Myz\n" {$Myx,$Myy,$Myz} >> source.txt
printf "%f\t%f\t%f\t\t\t\\\\\\ Mzx, Mzy, Mzz\n" {$Mzx,$Mzy,$Mzz} >> source.txt

printf "$nf\t$df\t$fim\t\t\\\\\\ frequency parameters\n" > fkparameters.txt
printf "$nu\t$du\t\t\t\\\\\\ slowness parameters\n" >> fkparameters.txt

## Compute wavelet and axes
expt --ofile expt.bin --nt $ntan --dt $dtan --a $omim --t0 $t0
rickerf --ofile waveletf.bin --nf $nf --df $df --t0 -$t0 --f0 $f0 --deriv 0 --ns 1 
linspace --ofile antaxis.bin --n $ntan --d  $dtan --o  $t0
linspace --ofile xaxis.bin --n $nrec --d $dxrec --o  $xrec0

################################################################################
## Run
################################################################################
fkm3d -o fk-f.bin -m medium.txt -f fkparameters.txt -s source.txt -r receivers.txt -v
freqmultfft --ofile fk-tr.bin --iarray fk-f.bin --ivector waveletf.bin --n2 $nf --df $df
prod --ofile fk-t.bin --ifile fk-tr.bin --ivector expt.bin

################################################################################
## Plot
################################################################################
figprfx="./Fig/"

aspectt=1.0
normalizet=0
tmin=-0.1
tmax=1.5

figname="fk"
comparewiggles --out "$figprfx$figname" --idata 'fk-t.bin ' --itaxes 'antaxis.bin' --ixaxis xaxis.bin --xlabel 'Time, seconds' --ylabel 'Field'  --normalize $normalizet --xmin $tmin --xmax $tmax --aspect $aspectt &