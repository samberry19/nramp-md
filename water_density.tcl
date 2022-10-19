#usage vmd -e water_density.tcl

package require volmap
package require pbctools
#load the psf and the dcd files and then center the trajectory
mol load psf system.psf
set list [lsort -dictionary [glob eq.*.dcd]]
set stp 1
foreach mol $list {
mol addfile $mol type dcd forst 0 last -1 step $stp filebonds 1 autobonds 1 waitfor all
}
pbc wrap -all -sel "not protein" -compound resid -center com -centersel "protein and name CA"


#align the protein to the first frame using the protein CA
set nf [molinfo top get numframes]
set all [atomselect top all]
set sel [atomselect top "protein and name CA"]
set ref [atomselect top "protein and name CA" frame 0]
for {set i 0} {$i<$nf} {incr i} {
$sel frame $i
$all frame $i
set fit [measure fit $sel $ref]
$all move $fit
}

#create the density of water molecules within 3.0 of protein. This writes a density file volmap_out.dx that can be loaded onto VMD. This density file is the occupancy file
#for the water molecules calculated for all the frames and then averaged. 

set wat [atomselect top "water and same residue as within 3.0 of protein"]

volmap occupancy $wat -res 1.0 -allframes -combine avg -o volmap_out.dx

quit
