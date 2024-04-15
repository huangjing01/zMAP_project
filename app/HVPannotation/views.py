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
from flask import render_template, request, jsonify,current_app,redirect, url_for,send_file,send_from_directory,flash
from . import HVPannotation
from werkzeug.utils import secure_filename

import os,string,random
import numpy as np
import pandas as pd
from jobs.HVPanno import HVPanno_
def processID(length=8,chars=string.ascii_letters+string.digits):
    return ''.join([random.choice(chars) for i in range(length)])


@HVPannotation.route('/', methods=['GET', 'POST'])
def index():
    return render_template('HVPannotate_index.html')



@HVPannotation.route('/result', methods=['POST'])
def result():
    condition = request.files['pvalue_results'] and request.files['z_statistic_matrix'] and request.files['sample_info']
    condition = condition and request.form.get('cluster_number_for_hypervariable') and request.form.get('minclustersize') and request.form.get('top') and request.form.get('cluster_number_for_top_proteins') and request.form.get('fdr')

    if condition:
        pvalue_res = request.files['pvalue_results']
        z_statistic = request.files['z_statistic_matrix']
        sample_info = request.files['sample_info']

        cluster_number_for_hypervariable = int(request.form.get('cluster_number_for_hypervariable'))
        minclustersize = int(request.form.get('minclustersize'))
        top = int(request.form.get('top'))
        cluster_number_for_top_proteins = int(request.form.get('cluster_number_for_top_proteins'))
        fdr = float(request.form.get('fdr'))

        proc_id = processID(8)
        os.chdir(current_app.config['UPLOAD_FOLDER'])
        while(os.path.exists(proc_id)):
            proc_id = processID(8)
        proc_id = proc_id + "_HVPanno"
        outdir = os.path.join(os.getcwd(),proc_id)
        os.mkdir(outdir)

        z_statistic_filename = secure_filename(z_statistic.filename)
        z_statistic.save(os.path.join(outdir,z_statistic_filename))

        pvalue_filename = secure_filename(pvalue_res.filename)
        pvalue_res.save(os.path.join(outdir,pvalue_filename))

        sample_info_filename = secure_filename(sample_info.filename)
        sample_info.save(os.path.join(outdir,sample_info_filename))



        task = HVPanno_.apply_async(kwargs={"pvalue_results" : os.path.join(outdir,pvalue_filename),
        "z_statistic_matrix" : os.path.join(outdir,z_statistic_filename),
        "sample_info" : os.path.join(outdir,sample_info_filename),
        "cluster_number_for_top_proteins" : cluster_number_for_top_proteins,
        "top" : top,
        "minclustersize" : minclustersize,
        "outdir" : outdir,
        "fdr" : fdr,
        "cluster_number_for_hypervariable": cluster_number_for_hypervariable})
            
        return render_template('HVPannotate_result.html',proc_id=proc_id,id=task.id)
    else:
        flash('please specify valid input files.')
        return redirect(url_for('HVPannotation.index'))


@HVPannotation.route('/status/<task_id>', methods=['POST','GET'])
def taskstatus(task_id):
    task = HVPanno_.AsyncResult(task_id)
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


@HVPannotation.route('/result/<filedir>', methods=['GET'])
def fetch_filedir(filedir):
    if filedir[-4:] == '.zip':
        return send_from_directory(os.path.abspath(current_app.config['UPLOAD_FOLDER']),filedir,as_attachment=True)
