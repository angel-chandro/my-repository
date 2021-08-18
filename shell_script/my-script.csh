#! /bin/tcsh -f
# nohup python3 [...].py >& out &

#set exec = 'cumu.py'
set exec = 'test.py'
#set exec = 'pklin.py'

set Testing = False

set nom = 'b'
set logpath = /home/$user/junk
set logname = ${logpath}/$nom.%A.%a.log
set job_file = ${logpath}/${nom}.job

if ($Testing == 'True') then
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
#SBATCH -A 32cores
#SBATCH -t 12:00:00
#SBATCH --mem=50000
#
echo $exec
python3 $exec
#
EOF

    sbatch $job_file
    rm $job_file
endif

echo 'All submited'
