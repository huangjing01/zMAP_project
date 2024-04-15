from flask import render_template, request, jsonify,current_app,redirect, url_for, send_from_directory,session,flash
from . import network
from werkzeug.utils import secure_filename

import os,string,random
import numpy as np
import pandas as pd
import json
from jobs.network import async_network
def processID(length=8,chars=string.ascii_letters+string.digits):
    return ''.join([random.choice(chars) for i in range(length)])

@network.route('/', methods=['GET', 'POST'])
def index():
    return render_template('network_index.html')

@network.route('/coexpression_network_example.xls', methods=['GET', 'POST'])
def example_data():
    return send_from_directory(os.path.abspath(os.path.join(current_app.config['UPLOAD_FOLDER'],'example')),'coexpression_network_example.xls',as_attachment=True)

@network.route('/validation', methods=['POST','GET'])
def validate(): 
#    if request.method == 'POST' and request.files['input-2']:
    if request.files['input-2']:
        file = request.files['input-2']
        if '.' in file.filename and file.filename.rsplit('.', 1)[1] == "xls":
            filename = secure_filename(file.filename)
            os.chdir(current_app.config['UPLOAD_FOLDER'])
            dir_name = processID(8)
            while(os.path.exists(dir_name)):
                dir_name = processID(8)
            os.mkdir(dir_name)
            file.save(os.path.join(os.getcwd(), dir_name,filename))
            return jsonify(error = "",id = dir_name)
        else:
            return jsonify(error = "Invalid filename suffix")
    else:
        return redirect(url_for("network.index"))



@network.route('/result/<process_id>', methods=['GET', 'POST'])
def result(process_id):
    if process_id[-4:] == '.zip':
        files_dir = os.path.join(current_app.config['UPLOAD_FOLDER'],process_id)
        return send_from_directory(os.path.abspath(current_app.config['UPLOAD_FOLDER']),process_id,as_attachment=True)
    else:
        os.chdir(os.path.join(current_app.config['UPLOAD_FOLDER'],process_id))
        enrichment = {}
        f2 = open("final_pcor_protein_label.txt","r")
        lines = f2.readlines()
        del(lines[0])
        for line in lines:
            label = line.split("\t")[1].replace("\n","")
            enrichment["clust"+str(label)]=[]
        f2.close()
        f1 = open("cluster_go_enrichment.txt","r")
        lines = f1.readlines()
        for line in lines:
            module = line.split("\t")[0]
            enrichment[module].append([line.split("\t")[1],line.split("\t")[2],line.split("\t")[3].replace("\n","")])
        f1.close()
        return render_template("network_result.html",enrichment=enrichment)

@network.route('/network_processing/<process_id>', methods=['POST','GET'])
def processing(process_id):
    files_dir = os.path.join(current_app.config['UPLOAD_FOLDER'],process_id)
    file = os.listdir(files_dir)[0]
    if request.values.get('method'):
        task = async_network.apply_async(kwargs={'filepath':os.path.join(files_dir,file),'module_detect':request.values.get('method')})
        return jsonify(error = "",process_id=process_id,method=request.values.get('method'),task_id= task.id)
    else:
        flash('please confirm parameters.')
        return redirect(url_for('network.index'))

@network.route('/status/<process_id>/<method>/<task_id>', methods=['POST','GET'])
def taskstatus(process_id,method,task_id):
    task = async_network.AsyncResult(task_id)
    response = {'state':'','status':[]}
    if task.state == 'PENDING':
        response['state'] = task.state
        response['status'].append('task Pending...')

    elif task.state != 'FAILURE': 
        response['state'] = task.state,
        response['status'] = task.info.get('status', '')

        if response['status'][-1] == 'Network finished':
            response['result'] = "over"
    else:
        response['state'] = task.state
        response['status'].append("Something Wrong.")
    return jsonify(response)

 
@network.route('/net_result/<process_id>', methods=['POST','GET'])
def net_result(process_id):
    all_json_res = []
    all_label = []
    os.chdir(os.path.join(current_app.config['UPLOAD_FOLDER'],process_id))
    f1 = open("final_pcor_protein_net.txt","r")
    lines = f1.readlines()
    del(lines[0])
    for line in lines:
        gene1 = line.split("\t")[0]
        gene2 = line.split("\t")[1].replace("\n","")
        edge = {"data":{"source":gene1,"target":gene2},"group":"edges"}
        all_json_res.append(edge)
    f1.close()
    f2 = open("final_pcor_protein_label.txt","r")
    lines = f2.readlines()
    del(lines[0])
    for line in lines:
        gene = line.split("\t")[0]
        label = line.split("\t")[1].replace("\n","")
        node = {"data":{"id":gene,"label":int(label),"parent":"clust"+str(label)},"group":"nodes"}
        all_json_res.append(node)
        all_label.append("clust"+str(label))
    for i in list(set(all_label)):
        node = {"data":{"id":i},"group":"nodes"}
        all_json_res.append(node)
    f2.close()
    return jsonify(all_json_res)
    

@network.route('/net_style/<process_id>', methods=['POST','GET'])
def net_style(process_id):
    os.chdir(os.path.join(current_app.config['UPLOAD_FOLDER'],process_id))
    style = [{"css":{"background-color":"rgb(137,208,245)","background-opacity":1.0,"border-color":"rgb(204,204,204)","border-opacity":1.0,"border-width":0.0,"color":"rgb(0,0,0)","content":"data(id)","font-family":"SansSerif.plain","font-size":12,"font-weight":"normal","height":35.0,"shape":"roundrectangle","text-halign":"center","text-opacity":1.0,"text-valign":"center","width":75.0},"selector":"node"},{"css":{"color":"rgb(0,0,0)","content":"","font-family":"Dialog.plain","font-size":10,"font-weight":"normal","line-color":"rgb(132,132,132)","line-style":"solid","opacity":1.0,"source-arrow-color":"rgb(0,0,0)","source-arrow-shape":"none","target-arrow-color":"rgb(0,0,0)","target-arrow-shape":"none","text-opacity":1.0,"width":2.0},"selector":"edge"},{"css":{"line-color":"rgb(255,0,0)"},"selector":"edge:selected"},{"css":{"background-opacity":"rgb(255,255,0)"},"selector":":parent"}]
    f = open("cluster_mapped_color.txt","r")
    lines = f.readlines()
    for line in lines:
        label = str(line.split("\t")[0])
        color = line.split("\t")[1].replace("\n","")
        print({"selector":"node[label = "+label+"]","css":{"background-color":color}})
        style.append({"selector":"node[label = "+label+"]","css":{"background-color":color}})
    f.close()
    style.append({"css":{"background-color":"rgb(255,255,0)"},"selector":"node:selected"})
    return jsonify(style)
