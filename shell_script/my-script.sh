#!/bin/sh
# nohup python3 [...].py >& out &

#set exec = 'cumu.py'
exec='test.py'
#set exec = 'pklin.py'
user=chandro

Testing=False

nom='archivo'
logpath=/home/$user/junk
logname=${logpath}/$nom.%A.%a.log
job_file=${logpath}/${nom}.job

if [ "$Testing" = True ]; then
    echo 'Testing'
    echo $exec
    python3 $exec
else
    cat > $job_file <<EOF
#! /bin/tcsh -ef
#
#SBATCH --ntasks 1 --cpus-per-task 16
#SBATCH -J ${nom}
#SBATCH -o ${logname}
#SBATCH -p all
#SBATCH -A 32cores
#SBATCH -t 12:00:00
#
echo $exec
python3 $exec
#
EOF

    sbatch $job_file
    rm $job_file
fi

echo 'All submited'
