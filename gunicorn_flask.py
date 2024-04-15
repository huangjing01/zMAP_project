import multiprocessing
 
bind = '127.0.0.1:9000'
#workers = multiprocessing.cpu_count() * 2 + 1
workers = 2
threads = 5
backlog = 2048
worker_class = "eventlet"
worker_connections = 1000
daemon = False
debug = True
proc_name = 'run'
pidfile = './logs/gunicorn.pid'
errorlog = './logs/gunicorn.log'

