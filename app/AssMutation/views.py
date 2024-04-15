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
from . import AssMutation
from werkzeug.utils import secure_filename

import os,string,random

from jobs.pQTL import async_pqtl
def processID(length=8,chars=string.ascii_letters+string.digits):
    return ''.join([random.choice(chars) for i in range(length)])


@AssMutation.route('/', methods=['GET', 'POST'])
def index():
    return render_template('AssMutation_index.html')


@AssMutation.route('/result', methods=['POST'])
def result():
    if request.files['z_statistic_matrix'] and request.files['mutation_f'] and request.files['covariates_f'] and request.files['gene_tss_location'] and request.files['chr_length']:

        mutation_f = request.files['mutation_f']
        z_statistic_matrix = request.files['z_statistic_matrix']
        covariates_f = request.files['covariates_f']
        gene_tss_location = request.files['gene_tss_location']
        chr_length = request.files['chr_length']
        fdr = float(request.form.get('fdr'))

        proc_id = processID(8)
        os.chdir(current_app.config['UPLOAD_FOLDER'])
        while(os.path.exists(proc_id)):
            proc_id = processID(8)
        proc_id = proc_id + "_pQTL"
        outdir = os.path.join(os.getcwd(),proc_id)
        os.mkdir(outdir)

        mutation_filename = secure_filename(mutation_f.filename)
        mutation_f.save(os.path.join(outdir,mutation_filename))

        z_statistic_filename = secure_filename(z_statistic_matrix.filename)
        z_statistic_matrix.save(os.path.join(outdir,z_statistic_filename))

        covariates_filename = secure_filename(covariates_f.filename)
        covariates_f.save(os.path.join(outdir,covariates_filename))

        genetss_filename = secure_filename(gene_tss_location.filename)
        gene_tss_location.save(os.path.join(outdir,genetss_filename))

        chr_filename = secure_filename(chr_length.filename)
        chr_length.save(os.path.join(outdir,chr_filename))


        task = async_pqtl.apply_async(kwargs={
        'mutation_f' : os.path.join(outdir,mutation_filename),
        'z_statistic_matrix' : os.path.join(outdir,z_statistic_filename),
        'covariates_f' : os.path.join(outdir,covariates_filename),
        'gene_tss_location' : os.path.join(outdir,genetss_filename),
        'chr_length' : os.path.join(outdir,chr_filename),
        'fdr' : fdr,
        'outdir' : outdir})
        
        return render_template('AssMutation_result.html',proc_id=proc_id,id=task.id)
    else:
        flash('please specify valid input files.')
        return redirect(url_for('AssMutation.index'))


@AssMutation.route('/status/<task_id>', methods=['POST','GET'])
def taskstatus(task_id):
    task = async_pqtl.AsyncResult(task_id)
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


@AssMutation.route('/result/<filedir>', methods=['GET'])
def fetch_filedir(filedir):
    if filedir[-4:] == '.zip':
        return send_from_directory(os.path.abspath(current_app.config['UPLOAD_FOLDER']),filedir,as_attachment=True)
