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

from flask import render_template, request, jsonify,current_app,redirect, url_for, send_from_directory,session,flash,send_file
from . import zMap
from jobs.zMap import async_zMap
import os,string,random,json,fnmatch
import numpy as np
import pandas as pd
import shutil

@zMap.route('/result/<process_id>', methods=['GET', 'POST'])
def result(process_id):  
    if process_id[-4:] == '.zip':
        files_dir = os.path.join(current_app.config['UPLOAD_FOLDER'],process_id)
        return send_from_directory(os.path.abspath(current_app.config['UPLOAD_FOLDER']),process_id,as_attachment=True)
    else:
        os.chdir(os.path.join(current_app.config['UPLOAD_FOLDER'],process_id))
        with open("zMAP_input.info.log","r") as g:
            a = json.load(g)
        g.close()
        task_id = a["task_id"]
        email = a["email"]
        return render_template("zMap_result.html",task_id=task_id,email=email)


@zMap.route('/image/<process_id>/<image_query>', methods=['GET', 'POST'])
def image_cls(process_id,image_query):
    os.chdir(os.path.join(current_app.config['UPLOAD_FOLDER'],process_id,'result'))
    image_list = []
    if image_query == "normalization":
        for path, dir, filelist in os.walk(os.getcwd()):
            for file_name in filelist:
                if file_name.endswith("_protein_intensity_boxplot.png"):
                    image_list.append(file_name)
        return jsonify(error="",id=image_list)
    elif image_query == "local":
        for path, dir, filelist in os.walk(os.getcwd()):
            for file_name in filelist:
                if file_name.endswith("_boxplot_for_r2_of_linear_regression.png"):
                    image_list.append(file_name)
        return jsonify(error="",id=image_list)
    elif image_query == "nonlinear_fitting":
        for path, dir, filelist in os.walk(os.getcwd()):
            for file_name in filelist:
                if file_name.endswith("model.png") and fnmatch.fnmatch(file_name,'*nonlinear*'):
                    image_list.append(file_name)
        return jsonify(error="",id=image_list)
    elif image_query == "theoretical_quantiles":
        for path, dir, filelist in os.walk(os.getcwd()):
            for file_name in filelist:
                if file_name.endswith("_theoretical_quantiles.png"):
                    image_list.append(file_name)
        return jsonify(error="",id=image_list)


@zMap.route('/result/<process_id>/<image>', methods=['GET', 'POST'])
def image_index(process_id,image):
    os.chdir(os.path.join(current_app.config['UPLOAD_FOLDER'],process_id,'result'))
    for path, dir, filelist in os.walk(os.getcwd()):
        if image in filelist:
            image_path = os.path.join(path,image)
            return send_file(image_path,mimetype='image/png')


def processID(length=8,chars=string.ascii_letters+string.digits):
    return ''.join([random.choice(chars) for i in range(length)])

@zMap.route('/', methods=['GET', 'POST'])
def index():
    return render_template('zMap_index.html')


@zMap.route('/validation/<input_type>', methods=['POST'])
def validate(input_type):
    try:
        with time_limit(2800):
            if request.method == 'POST':
                if input_type == 'protein_intensity':
                    file = request.files['input-1']
                    filename = "protein_intensity.txt"
                elif input_type == 'sample_info':
                    file = request.files['input-2']
                    filename = "sample_info.txt"
                else:
                    return jsonify(error="no data")
                os.chdir(current_app.config['UPLOAD_FOLDER'])
                id1 = processID(8)
                while(os.path.exists(id1+ '_' + filename)):
                    id1 = processID(8)

                file.save(os.path.join(os.getcwd(), id1+ '_' + filename))

                table_info = pd.read_csv(id1+ '_' + filename,index_col=0,sep="\t")

                if input_type == 'protein_intensity':
                    exp_gene_list = []
                    exp_sample_list = []
                    for index in table_info.index:
                        exp_gene_list.append(index)
                    if len(exp_gene_list) <= 1000:
                        return jsonify(error = "Not sufficent genes")
                    for column in table_info.columns:
                        if column.startswith('Unnamed'):
                            return jsonify(error = "Nans are included in the sample ID")
                        elif column in exp_sample_list:
                            return jsonify(error = "Conditions are not uniq")
                        else:
                            exp_sample_list.append(column)
                            if table_info[column].min() < 0:
                                return jsonify(error = "Negative value exists")
                    return jsonify(error = "", id = id1)   
                else:
                    if table_info.shape[0] > 0 and table_info.shape[1] > 0:
                        return jsonify(error = "", id = id1)
                    else:
                        return jsonify(error = "no sufficent sample information")
    except TimeoutException:
        return jsonify(error = "processing time exceeds limit")
    except:
        return jsonify(error = "please check your file format again")


@zMap.route('/zMap_processing/<process_id1>/<process_id2>', methods=['POST'])
def zMap_processing(process_id1,process_id2):
    intensity_file = os.path.join(current_app.config['UPLOAD_FOLDER'],process_id1 + '_protein_intensity.txt')
    sample_info_file = os.path.join(current_app.config['UPLOAD_FOLDER'],process_id2 + '_sample_info.txt')
    email = ''
    if request.form.get('window_size') and request.form.get('step_size') and request.form.get('percent') and request.form.get('method'):
        window_size = int(request.form.get('window_size'))
        step_size =int(request.form.get('step_size'))
        percent = int(request.form.get('percent'))/100
        method = request.form.get('method')
        if request.form.get('email'):
            email = request.form.get('email')
            session['email'] = email
        os.chdir(current_app.config['UPLOAD_FOLDER'])
        os.mkdir(process_id1)
        os.chdir(process_id1)
        intensity_file_path = os.path.join(current_app.config['UPLOAD_FOLDER'],process_id1,'protein_intensity.txt')
        sample_info_path = os.path.join(current_app.config['UPLOAD_FOLDER'],process_id1,'sample_info.txt')
        shutil.move(intensity_file, intensity_file_path)
        shutil.move(sample_info_file, sample_info_path)
        task = async_zMap.apply_async(kwargs={
            'intensity_file':intensity_file_path,
            'sample_info_file':sample_info_path,
            'window_size':window_size,
            'step_size':step_size,
            'percent':percent,
            'method':method,
            'email':email})
        
        zMap_input='zMAP_input.info.log'
        d = {'process_id':process_id1,'window_size':window_size,'step_size':step_size,'percent':percent,'method':method,'email':email,'task_id':task.id}
        with open(zMap_input,"w") as f:
            json.dump(d,f)
        f.close()
        return redirect(url_for("zMap.result",process_id=process_id1))
    else:
        flash('please confirm modeling parameters.')
        return redirect(url_for('zMap.index'))



@zMap.route('/status/<process_id>/<task_id>',methods=["POST","GET"])
def taskstatus(process_id,task_id):
    task = async_zMap.AsyncResult(task_id)
    response = {'state':'','status':[]}
    if task.state == 'PENDING':
        response['state'] = task.state
        response['status'].append('task pending...')

    elif task.state != 'FAILURE': 
        response['state'] = task.state,
        response['status'] = task.info.get('status', '')

        if response['status'][-1] == 'zMAP finished...':
            response['result'] = "over"
    else:
        # something went wrong in the background job
        response['state'] = task.state
        response['status'].append("Something Wrong.")
    return jsonify(response)


    

