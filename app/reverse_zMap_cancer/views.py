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
from . import reverse_zMap_cancer
from werkzeug.utils import secure_filename
from jobs.reverse_zMap_cancer import async_reverse_zMap_cancer
import os,string,random,json,fnmatch
import numpy as np
import pandas as pd


def processID(length=8,chars=string.ascii_letters+string.digits):
    return ''.join([random.choice(chars) for i in range(length)])



@reverse_zMap_cancer.route('/result/<process_id>', methods=['GET', 'POST'])
def result(process_id):  
    if process_id[-4:] == '.zip':
        files_dir = os.path.join(current_app.config['UPLOAD_FOLDER'],process_id)
        return send_from_directory(os.path.abspath(current_app.config['UPLOAD_FOLDER']),process_id,as_attachment=True)
    else:
        os.chdir(os.path.join(current_app.config['UPLOAD_FOLDER'],process_id))
        with open("reverse_zMAP_cancer_input.info.log","r") as g:
            a = json.load(g)
        g.close()
        task_id = a["task_id"]
        email = a["email"]
        return render_template("reverse_zMap_cancer_result.html",task_id=task_id,email=email)

@reverse_zMap_cancer.route('/image/<process_id>/<image_query>', methods=['GET', 'POST'])
def image_cls(process_id,image_query):
    os.chdir(os.path.join(current_app.config['UPLOAD_FOLDER'],process_id,'result'))
    image_list = []
    if image_query == "qc":
        for path, dir, filelist in os.walk(os.getcwd()):
            for file_name in filelist:
                if file_name in ["qc_table.png","qc_boxplot.png"]:
                    image_list.append(file_name)
        return jsonify(error="",id=image_list)
    elif image_query == "distribution":
        for path, dir, filelist in os.walk(os.getcwd()):
            for file_name in filelist:
                if file_name in ["clustermap.png","PC1_PC2_scatterplot.png","distribution_of_m_z.png","cumulation.png"]:
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


@reverse_zMap_cancer.route('/result/<process_id>/<image>', methods=['GET', 'POST'])
def image_index(process_id,image):
    os.chdir(os.path.join(current_app.config['UPLOAD_FOLDER'],process_id,'result'))
    for path, dir, filelist in os.walk(os.getcwd()):
        if image in filelist:
            image_path = os.path.join(path,image)
            return send_file(image_path,mimetype='image/png')






@reverse_zMap_cancer.route('/reverse_zMAP.cancer_example_data.xlsx', methods=['GET', 'POST'])
def example_data():
    return send_from_directory(os.path.abspath(os.path.join(current_app.config['UPLOAD_FOLDER'],'example')),'reverse_zMAP.cancer_example_data.xlsx',as_attachment=True)




@reverse_zMap_cancer.route('/', methods=['GET', 'POST'])
def index():
    return render_template('reverse_zMap_cancer_index.html')

#file content should be validated

@reverse_zMap_cancer.route('/validation', methods=['POST'])
def validate():
    try:
        with time_limit(2800):
            if request.method == 'POST' and request.files['input-1']: 
                file = request.files['input-1']
                if '.' in file.filename and file.filename.rsplit('.', 1)[1] == "xlsx":
                    filename = secure_filename(file.filename)
                    os.chdir(current_app.config['UPLOAD_FOLDER'])
                    dir_name = processID(8)
                    while(os.path.exists(dir_name)):
                        dir_name = processID(8)
                    os.mkdir(dir_name)
                    os.chdir(dir_name)
                    file.save(os.path.join(os.getcwd(), filename))
                    raw_intensity_data = pd.ExcelFile(filename)
                    sheet_names = raw_intensity_data.sheet_names
                    intensity_df = raw_intensity_data.parse(sheet_names[0],index_col=0)
                    batch_df = raw_intensity_data.parse(sheet_names[1],index_col=0)

                    exp_gene_list = []
                    exp_sample_list = []
                    g_sample_list = []
                    for index in intensity_df.index:
                        if isinstance(index,float):
                            return jsonify(error = "Nans are included in the gene ID")
                        elif index in exp_gene_list:
                            return jsonify(error = "Gene IDs were not uniq")
                        else:
                            exp_gene_list.append(index)
                    if len(exp_gene_list) <= 1000:
                        return jsonify(error = "Not sufficent genes")
                    for column in intensity_df.columns:
                        if column.startswith('Unnamed'):
                            return jsonify(error = "Nans are included in the sample ID")
                        elif column in exp_sample_list:
                            return jsonify(error = "Sample IDs were not uniq")
                        elif (column.startswith("T") or column.startswith("N")) or column.endswith("MIX"):
                            exp_sample_list.append(column)
                            if intensity_df[column].isnull().sum() > 0:
                                return jsonify(error = "Nans are included in expression values")
                            elif intensity_df[column].min() < 0:
                                return jsonify(error = "Negative values exists")
                            else:
                                pass
                        else:
                            return jsonify(error = "Sample id should starts with T(tumor sample) or N(normal sample),mix sample id should ends with(MIX)")
                    for column in batch_df.columns:
                        g_sample_list.extend(batch_df[column][batch_df[column].notnull()].tolist())
                    if set(exp_sample_list) != set(g_sample_list):
                        return jsonify(error = "Sample list not exactly the same")
                else:
                    return jsonify(error = "Invalid filename suffix")
                return jsonify(error = "",id = dir_name)  
            else:
                return redirect(url_for("reverse_zMap_cancer.index"))
    except TimeoutException:
        return jsonify(error = "processing time exceeds limit")
    except:
        return jsonify(error = "please check your file format again")            



@reverse_zMap_cancer.route('/reverse_zMap_cancer_processing/<process_id>', methods=['POST'])
def reverse_zMap_cancer_processing(process_id):
    files_dir = os.path.join(current_app.config['UPLOAD_FOLDER'],process_id)
    file = os.listdir(files_dir)[0]
    email = ''
    if request.form.get('window_size') and request.form.get('step_size') and request.form.get('percent') and request.form.get('method'):
        window_size = int(request.form.get('window_size'))
        step_size =int(request.form.get('step_size'))
        percent = int(request.form.get('percent'))/100
        method = request.form.get('method')
        if request.form.get('email'):
            email = request.form.get('email')
            session['email'] = email
        task = async_reverse_zMap_cancer.apply_async(kwargs={'filepath':os.path.join(files_dir,file),'window_size':window_size,'step_size':step_size,'percent':percent,'method':method,'email':email})
        os.chdir(files_dir)

        reverse_zMap_cancer_input='reverse_zMAP_cancer_input.info.log'
        d = {'filepath':os.path.join(files_dir,file),'process_id':process_id,'window_size':window_size,'step_size':step_size,'percent':percent,'method':method,'email':email,'task_id':task.id}
        with open(reverse_zMap_cancer_input,"w") as f:
            json.dump(d,f)
        f.close()
        return redirect(url_for("reverse_zMap_cancer.result",process_id=process_id))
    else:
        flash('please confirm modeling parameters.')
        return redirect(url_for('reverse_zMap_cancer.index'))


@reverse_zMap_cancer.route('/status/<process_id>/<task_id>',methods=["POST","GET"])
def taskstatus(process_id,task_id):
    task = async_reverse_zMap_cancer.AsyncResult(task_id)
    response = {'state':'','status':[]}
    if task.state == 'PENDING':
        response['state'] = task.state
        response['status'].append('task Pending...')

    elif task.state != 'FAILURE': 
        response['state'] = task.state,
        response['status'] = task.info.get('status', '')

        if response['status'][-1] == 'reverse_zMAP.cancer finished':
            response['result'] = "over"
    else:
        # something went wrong in the background job
        response['state'] = task.state
        response['status'].append("Something Wrong.")
    return jsonify(response)






