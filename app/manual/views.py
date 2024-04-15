from flask import render_template, request,current_app,redirect,url_for,send_file,send_from_directory
from . import manual
import os


@manual.route('/', methods=['GET', 'POST'])
def home():
    return render_template('zMAP_paper.html')

@manual.route('/example/<example_file>', methods=['GET', 'POST'])
def example_data(example_file):
    if example_file in os.listdir(os.path.join(current_app.config['UPLOAD_FOLDER'],'example/input')):
        if example_file.endswith('.txt'):
            if example_file.find('sample_info') != -1:
                return send_from_directory(os.path.abspath(os.path.join(current_app.config['UPLOAD_FOLDER'],'example/input')),example_file,as_attachment=False)
            else:
                return send_from_directory(os.path.abspath(os.path.join(current_app.config['UPLOAD_FOLDER'],'example/input')),example_file,as_attachment=True)
        elif example_file.endswith(".pdf"):
            return send_file(os.path.join(os.path.abspath(os.path.join(current_app.config['UPLOAD_FOLDER'],'example/input')),example_file), mimetype='application/pdf')
        else:
            flash("No example file you specified")
    else:
        flash("Your queried file doesn't exist")
    

@manual.route('/example/<image>', methods=['GET', 'POST'])
def image_index(process_id,image):
    os.chdir(os.path.join(current_app.config['UPLOAD_FOLDER'],process_id,'result'))
    for path, dir, filelist in os.walk(os.getcwd()):
        if image in filelist:
            image_path = os.path.join(path,image)
            return send_file(image_path,mimetype='image/png')

