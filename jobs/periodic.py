# -*- coding: utf-8 -*-
"""
async task: add_periodic_task
"""

from app import celery
from celery.schedules import crontab
from flask import current_app,Flask
import os,time,datetime,re,shutil
user_upload_dir = os.path.abspath(os.path.dirname(__file__))
@celery.task(name='periodic')
def test():
    os.chdir(os.path.join(os.path.dirname(os.path.dirname(__file__)),"uploads"))
    curr_time = time.localtime()
    curr_time_day = datetime.datetime(curr_time.tm_year,curr_time.tm_mon,curr_time.tm_mday)
    if(len(os.listdir(os.getcwd())) != 1):
        for f in os.listdir(os.getcwd()):
            if f !="example":
                dir_c_time = time.localtime(os.path.getmtime(f))
                dir_time_day =  datetime.datetime(dir_c_time.tm_year,dir_c_time.tm_mon,dir_c_time.tm_mday)
                exist_days = int(re.split('[, :]',str(curr_time_day - dir_time_day))[0])
                if exist_days > 0:
                    filepath = os.path.join(os.getcwd(),f)
                    if os.path.isfile(filepath):
                        os.remove(filepath)
                    elif os.path.isdir(filepath):
                        shutil.rmtree(filepath)




		
@celery.on_after_finalize.connect
def setup_periodic_tasks(sender, **kwargs):
    # Calls test('hello') every 10 seconds.
    sender.add_periodic_task(60*60*24*7, test.s(), name='periodic')
