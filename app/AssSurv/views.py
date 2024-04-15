from __future__ import with_statement
import signal
from contextlib import contextmanager

class TimeoutException(Exception): pass

@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
from flask import render_template, request, jsonify,current_app,redirect, url_for,send_from_directory,flash
from . import AssSurv
from werkzeug.utils import secure_filename

import os,string,random
import numpy as np
import pandas as pd
from jobs.AssSurv import cox_regression


def processID(length=8,chars=string.ascii_letters+string.digits):
    return ''.join([random.choice(chars) for i in range(length)])


@AssSurv.route('/', methods=['GET', 'POST'])
def index():
    return render_template('AssSurv_index.html')


@AssSurv.route('/result', methods=['POST'])
def result():
    if request.files['input-1'] and request.files['input-2']:

        z_statistic_f = request.files['input-1']
        survival_df = request.files['input-2']

        proc_id = processID(8)
        os.chdir(current_app.config['UPLOAD_FOLDER'])
        while(os.path.exists(proc_id)):
            proc_id = processID(8)
        proc_id = proc_id + "_AssSurv"
        outdir = os.path.join(os.getcwd(),proc_id)
        os.mkdir(outdir)

        z_statistic_filename = secure_filename(z_statistic_f.filename)
        survival_df_filename = secure_filename(survival_df.filename)
        z_statistic_f.save(os.path.join(outdir,z_statistic_filename))
        survival_df.save(os.path.join(outdir,survival_df_filename))

        task = cox_regression.apply_async(kwargs={'z_statistic':os.path.join(outdir,z_statistic_filename),'survival_f':os.path.join(outdir,survival_df_filename), 'outdir': outdir})

        return render_template('AssSurv_result.html',proc_id=proc_id,id=task.id)
    else:
        flash('please specify valid input files.')
        return redirect(url_for('AssSurv.index'))

@AssSurv.route('/status/<task_id>', methods=['POST','GET'])
def taskstatus(task_id):
    task = cox_regression.AsyncResult(task_id)
    response = {'state':'','status':[]}
    if task.state == 'PENDING':
        response['state'] = task.state
        response['status'].append('task Pending...')

    elif task.state != 'FAILURE': 
        response['state'] = task.state,
        response['status'] = task.info.get('status', '')

        if response['status'][-1] == 'All finished.':
            response['result'] = "over"
    else:
        response['state'] = task.state
        response['status'].append("Something Wrong.")
    return jsonify(response)


@AssSurv.route('/result/<filedir>', methods=['GET'])
def fetch_filedir(filedir):
    if filedir[-4:] == '.zip':
        return send_from_directory(os.path.abspath(current_app.config['UPLOAD_FOLDER']),filedir,as_attachment=True)
