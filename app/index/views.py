from flask import render_template, request,current_app,redirect,url_for
from . import index
import os


@index.route('/', methods=['GET', 'POST'])
def home():
    return render_template('index.html')

@index.route('/query', methods=['GET', 'POST'])
def query():
    process_id = request.form.get('query_id')
    os.chdir(current_app.config['UPLOAD_FOLDER'])
    if os.path.exists(process_id):
        os.chdir(process_id)
        if os.path.isfile("zMAP_input.info.log"):
            return redirect(url_for("zMap.result",process_id=process_id))
        elif os.path.isfile("reverse_zMap_input.info.log"):
            return redirect(url_for("reverse_zMAP.result",process_id=process_id))
        else:
            pass
    else:
        return "The process id you entered doesn't exists."
