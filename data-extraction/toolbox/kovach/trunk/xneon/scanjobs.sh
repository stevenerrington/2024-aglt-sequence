#!/bin/bash
#
# A script to scan for new matlab jobs submitted through xneon.
#
# Use:
#
#   ./scanjobs.sh [directory]
# 
# Check directory (default ~/jobs) for new jobs every 10 second.
#
#

# Christopher Kovach 2016


## This is the directory that will be scanned.
jobsdir=${1:-~/xneon_submit/jobs}

## Repeat after this many seconds
repeat_every=10

## Delete jobs after they have been idle (created but not submitted) for this long.
idlelimit=3600


cd $jobsdir


[ -e "${jobsdir}/scanned" ] && scanfi=$(cat ${jobsdir}/scanned)

isscanning=$( echo $scanfi | awk '{print $1}')
scanpid=$( echo $scanfi | awk '{print $2}')
#echo $isscanning
## Don't run if another scan is still going
if [[ ${isscanning:-done} != done ]]
then
    rm ${jobsdir}/scanned
    echo Aborting...
    exit 101
elif  [ "$isscanning" ] && [ $(($(date +%s)-$(date -r ${jobsdir}/scanned +%s ))) -lt $((2*$repeat_every)) ] && [[ ${scanpid:-0} -ne $$ ]] 
then
    echo $0 appears to be running already or tried to run and crashed
    exit 100
fi



echo "running $$" > ${jobsdir}/scanned

diskuse=$(df  ~ | grep -E -o '[0-9]+%')
echo "disk use: $diskuse" >> ${jobsdir}/scanned



jobs=$( ls -d ${jobsdir}/*/ 2>/dev/null )

for job in ${jobs}
do
    
    jobdir=${job}
    chmod 777 ${jobdir}
    statusfile=${jobdir}/status/current

    if ! [ -e $statusfile ] 
    then
        echo skipping ${jobdir}
        continue
    fi

    state=$(cat $statusfile )

    [ -e ${jobdir}/status/jobid ] && jobid=$(cat ${jobdir}/status/jobid)  
 
    case "${state}" in
        finished | error)
            continue
            ;;

        jobover)
            qacct -j ${jobid} > ${jobdir}/status/qacct
            echo finished > ${statusfile}
            ;;

        submitted | running | waiting | suspended | error | deleting)
            
            # Get the status information and parse the result
            qst=$(qstat | grep ${jobid} 2>&1 )
            [ "$qst" ] || echo finished > ${statusfile} #&& qacct -j ${jobid} > ${jobdir}/status/qacct
            
            echo "${qst}" > ${jobdir}/status/qstat
            case "$(echo ${qst} | awk '{print $5}')" in
                *[r]* )
                    echo running > ${statusfile}
                    ;;
                *[w]*)
                    echo waiting > ${statusfile}
                    ;;
                *[sS]*) echo suspended > ${statusfile}
                    ;;
                *[E]*) echo error > ${statusfile}
                    ;;
                *[d]*) echo deleting > ${statusfile}
                    ;;
             esac
             ;;
        closed)
            [ -e "${jobdir}" ] && rm -rf ${jobdir}
             qdel $jobid 
            ;;

        open)
            
            ## Submit an open request to the queue through qsub
            qsb=$(qsub -o ${jobdir}/log/output -e ${jobdir}/log/error ${jobdir}/bash/submit.sh)

            ## Check if there was an error and update status
            ( (( $? )) ||  [ $( echo submitted > ${statusfile} ) ] ) && echo error > ${statusfile}

            ## get the job ID number from the output
            echo $qsb | sed 's/Your job[-ary]* \([0-9]*\).*/\1/' > ${jobdir}/status/jobid
            
            qst=$(qstat | grep ${jobid} 2>&1 )
            echo "${qst}" > ${jobdir}/status/qstat
            ;;
        *)

           ## In all other cases delete the job if there has been no change
           ## of status for an hour or more.
           [ $(($(date +%s)-$(date -r $statusfile +%s ))) -lt $idlelimit ] || rm -rf ${jobdir}
    esac
done

echo "done $$" > ${jobsdir}/scanned

diskuse=$(df  ~ | grep -E -o '[0-9]+%')
echo "disk use: $diskuse" >> ${jobsdir}/scanned

### Rotate the log file
maxlogsize=2M
keepNlogs=10
printf "\ncompress \n\n\"./jobs.log\" {\n\trotate $keepNlogs\\n\tsize $maxlogsize\n\t\n}\n" > .logr.conf
[ -e jobs.log ] && /usr/sbin/logrotate -s .logr.status .logr.conf

### Sleep for repeat_every seconds and repeat
sleep $repeat_every && exec $0
