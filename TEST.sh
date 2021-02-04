kill -9 `ps -ef|grep ML|awk '{print $2}'`
rm nohup.out
nohup bash RUN_ML.sh &
