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
from . import AssClinical
from werkzeug.utils import secure_filename

import os,string,random

from jobs.AssClinical import association_with_clinical_feature_
def processID(length=8,chars=string.ascii_letters+string.digits):
    return ''.join([random.choice(chars) for i in range(length)])


@AssClinical.route('/', methods=['GET', 'POST'])
def index():
    return render_template('AssClinical_index.html')


@AssClinical.route('/result', methods=['POST'])
def result():
    condition = request.files['input-1'] and request.files['input-2'] and request.files['input-3'] and request.files['input-4'] and request.files['input-5']
    condition = condition and request.form.get('discrete') and request.form.get('continuous') and request.form.get('fdr') and request.form.get('cluster_n')

    if condition:
        z_statistic = request.files['input-1']
        cluster_f = request.files['input-2']
        clinical_info = request.files['input-3']
        color_f = request.files['input-4']
        colorbar_f = request.files['input-5']

        discrete = request.form.get('discrete')
        continuous = request.form.get('continuous')
        fdr = float(request.form.get('fdr'))
        cluster_n = int(request.form.get('cluster_n'))

        proc_id = processID(8)
        os.chdir(current_app.config['UPLOAD_FOLDER'])
        while(os.path.exists(proc_id)):
            proc_id = processID(8)
        proc_id = proc_id + "_AssClincal"
        outdir = os.path.join(os.getcwd(),proc_id)
        os.mkdir(outdir)

        z_statistic_filename = secure_filename(z_statistic.filename)
        z_statistic.save(os.path.join(outdir,z_statistic_filename))

        cluster_filename = secure_filename(cluster_f.filename)
        cluster_f.save(os.path.join(outdir,cluster_filename))

        clinical_filename = secure_filename(clinical_info.filename)
        clinical_info.save(os.path.join(outdir,clinical_filename))

        colorlist_filename = secure_filename(color_f.filename)
        color_f.save(os.path.join(outdir,colorlist_filename))

        colobar_filename = secure_filename(colorbar_f.filename)
        colorbar_f.save(os.path.join(outdir,colobar_filename))


        task = association_with_clinical_feature_.apply_async(kwargs={"z_statistic_matrix" : os.path.join(outdir,z_statistic_filename),
        "cluster" : os.path.join(outdir,cluster_filename),
        "clinical_info" : os.path.join(outdir,clinical_filename),
        "discrete" : discrete.strip().split(','),
        "continuous" : continuous.strip().split(','),
        "color_f" : os.path.join(outdir,colorlist_filename),
        "colorbar_f" : os.path.join(outdir,colobar_filename),
        "outdir" : outdir,
        "fdr" : fdr,
        "cluster_n": cluster_n})
            
        return render_template('AssClinical_result.html',proc_id=proc_id,id=task.id)
    else:
        flash('please specify valid input files.')
        return redirect(url_for('AssClinical.index'))


@AssClinical.route('/status/<task_id>', methods=['POST','GET'])
def taskstatus(task_id):
    task = association_with_clinical_feature_.AsyncResult(task_id)
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


@AssClinical.route('/result/<filedir>', methods=['GET'])
def fetch_filedir(filedir):
    if filedir[-4:] == '.zip':
        return send_from_directory(os.path.abspath(current_app.config['UPLOAD_FOLDER']),filedir,as_attachment=True)
