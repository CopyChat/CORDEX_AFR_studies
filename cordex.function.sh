#!/bin/bash - 
#===============================================================================
#
#          FILE: function.sh
# 
#         USAGE: ./function.sh 
# 
#   DESCRIPTION:  
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Tang (Tang), tangchao90908@sina.com
#  ORGANIZATION: KLA
#       CREATED: 02/02/17 16:41:10 RET
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

function copydaily()
{
	ln -sf /worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/daily/rsds/*/nc ./
}

function copymonthly()
{
	#echo $1
	ln -sf /worktmp2/RegCM_DATA/MODEL_DATA/CORDEX/AFR-44/daily/rsds/*/*hist-rcp85*nc ./
}

function selyear()
{
	for f in *calendar.nc
	do
		echo $f
		echo $(eval seq -s "," $1 $2)
		cdo -b 64 selyear,$(eval seq -s "," $1 $2) $f $f.${1}-${2}.nc
		cdo -b 64 selyear,$(eval seq -s "," $3 $4) $f $f.${3}-${4}.nc
	done
}

function selregion()
{
	for f in *calendar.nc.????-????.nc
	do
		echo $f

		cdo -b 64 sellonlatbox,10,30,0,-10   $f $f.1.nc
		cdo -b 64 sellonlatbox,30,40,0,-10   $f $f.2.nc
		cdo -b 64 sellonlatbox,15,40,-10,-15 $f $f.3.nc
		cdo -b 64 sellonlatbox,15,20,-15,-30 $f $f.4.nc
		cdo -b 64 sellonlatbox,20,35,-15,-35 $f $f.5.nc
		cdo -b 64 sellonlatbox,45,50,-15,-25 $f $f.6.nc

		cdo -b 64 sellonlatbox,54.9,58,-19.5,-21.5 $f $f.7.nc # RUN & MAU

    	cdo -b 64 sellonlatbox,54.9,56,-20.5,-21.5 $f $f.8.nc # RUN

    	cdo -b 64 sellonlatbox,57.2,58,-19.5,-20.8 $f $f.9.nc # MAU

	echo "move to next one ..."
	done
}

function fldmean()
{
	for f in *nc.?.nc
	do
		echo $f
		cdo -b 64 fldmean $f $f.fldmean.nc
	done
}

function detrending()
{
	for f in *fldmean.nc
	do
		echo $f

		cdo detrend $f $f.detrend.nc

		echo "move to next one ..."
	done
}
function maskannual()
{
	for f in *detrend.nc
	do

		#cdo -b 64 ydaymean $f $f.ydaymean.nc
		#cdo -b 64 ydaysub $f $f.ydaymean.nc $f.daymean.maskannual.nc

		#cdo -b 64 monmean $f $f.monmean.nc
		#cdo -b 64 ymonmean $f $f.ymonmean.nc
		#cdo -b 64 ymonsub $f.monmean.nc  $f.ymonmean.nc $f.monmean.maskannual.nc

		cdo -b 64 yearmean $f $f.yearmean.masknon.nc
	done
}
#function maskannual2()
#{
#for i in 7 8 9
#do
	#for f in *nc.$i.nc.fldmean.nc.detrend.nc
	#do

		#echo $f

		#cdo -b 64 ydaymean $f $f.ydaymean.nc
		#cdo -b 64 ydaysub $f $f.ydaymean.nc $f.daymean.maskannual.nc

		#cdo -b 64 monmean $f $f.monmean.nc
		#cdo -b 64 ymonmean $f $f.ymonmean.nc
		#cdo -b 64 ymonsub $f.monmean.nc  $f.ymonmean.nc $f.monmean.maskannual.nc

		#cdo -b 64 yearmean $f $f.yearmean.masknon.nc
	#done
#done
#}
#=================================================== cal
function minmax()
{
for i in 1 2 3 4 5 6 7 8 9
do
	for f in $(cat list)
	do

		echo $i $f
		cdo sub $f.calendar.nc.$3-$4.nc.$i.nc $f.calendar.nc.$1-$2.nc.$i.nc $f.calendar.daychanges.$i.nc

		cdo timmin $f.calendar.daychanges.$i.nc $f.calendar.daychanges.timmin.$i.nc 
		cdo timmax $f.calendar.daychanges.$i.nc $f.calendar.daychanges.timmax.$i.nc 

		#rsds_AFR-44_NCC-NorESM1-M_SMHI-RCA4_v1.hist_rcp85.day.1951-2100.nc.calendar.nc.2041-2070.nc.6.nc
	done

	#rsds_AFR-44_CCCma-CanESM2_SMHI-RCA4_v1.hist_rcp85.day.1951-2100.nc.calendar.daychanges.timmax.8.nc
	cdo ensmin *calendar.daychanges.timmin.$i.nc rsds_AFR-44.calendar.daychanges.timmin.ensmin.$i.nc
	cdo ensmax *calendar.daychanges.timmax.$i.nc rsds_AFR-44.calendar.daychanges.timmax.ensmax.$i.nc
done

}
function meanchanges()
{

for i in 1 2 3 4 5 6 7 8 9
do
    #for f in $(ls *1976*fldmean.nc); do echo ${f%.1976*}; done | uniq > original.cal

    echo -------- $i -----------

    for f in $(cat original.cal)
    do
        cdo sub $f.$3-$4.nc.$i.nc.fldmean.nc $f.$1-$2.nc.$i.nc.fldmean.nc $f.fldmean.daychanges.$i.nc
        cdo timmean $f.fldmean.daychanges.$i.nc $f.fldmean.daychanges.timmean.$i.nc
    done
    cdo ensmean *.fldmean.daychanges.timmean.$i.nc rsds_AFR-44.fldmean.daychanges.timmean.ensmean.$i.nc
done
}


function calstd()
{

    for f in *fldmean.nc.detrend.nc.*mask*
    do
        echo $f 
        cdo timstd $f $f.timstd.nc
    done
}
function meanstd()
{
    for i in 1 2 3 4 5 6 7 8 9
    do
        echo 000000000000000 $i 00000000000000
        for t in daymean monmean yearmean
        do
            echo 000000000000000 $t 00000000000000
            echo *nc.2041-2070.nc.$i.nc.fldmean.nc.*timstd.nc
            cdo ensmean *nc.$1-$2.nc.$i.nc.*$t.*.timstd.nc rsds_AFR-44.$1-$2.$i.$t.timstd.ensstd.nc
            cdo ensmean *nc.$3-$4.nc.$i.nc.*$t.*.timstd.nc rsds_AFR-44.$3-$4.$i.$t.timstd.ensstd.nc
        done
    done
}

function printstd()
{
    # for yearmean
    for i in 1 2 3 4 5 6 7 8 9
    do
        echo $i
        for t in yearmean
        do
            echo $t
            for f in *.nc.$1-$2.nc.$i.nc.*.$t.masknon.nc.timstd.nc
            do
                refstd=$(ncks -v rsds $f | grep "W m-2" | grep time | awk -F "=" '{printf "%5.4f\n", $5}')
                
                futurefile=${f%.nc.$1-$2*}.nc.$3-$4.nc.$i.nc.*.$t.masknon.nc.timstd.nc
                #echo $futurefile
                for k in $(ls $futurefile)
                do
                    fustd=$(ncks -v rsds $k | grep "W m-2" | grep time | awk -F "=" '{printf "%5.4f\n", $5}')
                done

                echo $refstd $fustd | awk '{printf "%5.3f\n",$2/$1}'
            done
        done
        
        # for ensstd
        refstd2=$(ncks -v rsds rsds_AFR-44.$1-$2.$i.$t.timstd.ensstd.nc | grep "W m-2" | grep time | awk -F "=" '{printf "%5.4f\n", $5}')

        fustd2=$(ncks -v rsds rsds_AFR-44.$3-$4.$i.$t.timstd.ensstd.nc | grep "W m-2" | grep time | awk -F "=" '{printf "%5.4f\n", $5}')

        echo $refstd2 $fustd2| awk '{printf "%5.3f\n",$2/$1}'
    done
#=================================================== 
    # for daily and monthly
    for i in 1 2 3 4 5 6 7 8 9
    do
        echo $i
        for t in daymean monmean
        do
            echo $t
            for f in *.nc.$1-$2.nc.$i.nc.*.$t.maskannual.nc.timstd.nc
            do
                refstd=$(ncks -v rsds $f | grep "W m-2" | grep time | awk -F "=" '{printf "%5.4f\n", $5}')
                
                futurefile=${f%.nc.$1-$2*}.nc.$3-$4.nc.$i.nc.*.$t.maskannual.nc.timstd.nc
                #echo $futurefile
                for k in $(ls $futurefile)
                do
                    fustd=$(ncks -v rsds $k | grep "W m-2" | grep time | awk -F "=" '{printf "%5.4f\n", $5}')
                done

                echo $refstd $fustd | awk '{printf "%5.3f\n",$2/$1}'
            done
        done
        
        # for ensstd
        refstd2=$(ncks -v rsds rsds_AFR-44.$3-$4.$i.$t.timstd.ensstd.nc | grep "W m-2" | grep time | awk -F "=" '{printf "%5.4f\n", $5}')

        fustd2=$(ncks -v rsds rsds_AFR-44.$1-$2.$i.$t.timstd.ensstd.nc | grep "W m-2" | grep time | awk -F "=" '{printf "%5.4f\n", $5}')

        echo $refstd2 $fustd2 | awk '{printf "%5.3f\n",$2/$1}'
    done
}
function mean_ref()
{
for i in 1 2 3 4 5 6 7 8 9
do
    #for f in *day.1951-2100.nc.calendar.nc.1976-2005.nc.$i.nc.fldmean.nc
    #do
        #echo $f
        #cdo timmean $f $f.timmean.nc
    #done
    cdo ensmean *day.1951-2100.nc.calendar.nc.1976-2005.nc.$i.nc.fldmean.nc.timmean.nc rsds_AFR.1976-2005.fldmean.timmean.ensmean.$i.nc
done

}
function print_mean_ref()
{
for i in 1 2 3 4 5 6 7 8 9
do
    echo -------- $i -----------
    for f in *day.1951-2100.nc.calendar.nc.1976-2005.nc.$i.nc.fldmean.nc.timmean.nc
    do
        meanref=$(ncks -v rsds $f | grep "W m-2" | grep time | awk -F "=" '{printf "%5.4f\n", $5}')
        echo $meanref
    done
done
}
function mean_future()
{
    for f in *day.1951-2100.nc.calendar.nc.2041-2070.nc.*.nc.fldmean.nc
    do
        echo $f
        cdo timmean $f $f.timmean.nc
    done

}
function print_mean_future()
{
for i in 1 2 3 4 5 6 7 8 9
do
    echo -------- $i -----------
    for f in *day.1951-2100.nc.calendar.nc.2041-2070.nc.$i.nc.fldmean.nc.timmean.nc
    do
        meanref=$(ncks -v rsds $f | grep "W m-2" | grep time | awk -F "=" '{printf "%5.4f\n", $5}')
        echo $meanref
    done
done
}

function printmeanchanges()
{

    #           region_1 region_2 ... region_i ...
    #  model_1
    #  model_2
    #  ...
    #  model_m
    #  ...

    for m in $(cat original.cal)
    do
        for i in 1 2 3 4 5 6 7 8 9
        do
            meanchanges=$(ncks -v rsds $m.fldmean.daychanges.timmean.$i.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            mean_ref=$(ncks -v rsds $m.1976-2005.nc.$i.nc.fldmean.nc.timmean.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            echo $meanchanges $mean_ref | awk '{printf "%5.3f\t",$1/$4}'
        done
        echo ""
    done

        for i in 1 2 3 4 5 6 7 8 9
        do
            ensmeanchanges=$(ncks -v rsds rsds_AFR-44.fldmean.daychanges.timmean.ensmean.$i.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            ensmean_ref=$(ncks -v rsds rsds_AFR.1976-2005.fldmean.timmean.ensmean.$i.nc | grep "W m-2" | grep time | awk -F "=" '{print $5}')
            echo $ensmeanchanges $ensmean_ref | awk '{printf "%5.3f\t",$1/$4}'
        done
        echo ""
}
