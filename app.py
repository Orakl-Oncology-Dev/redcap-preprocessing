from flask import Flask, render_template, request, send_file
import os
from werkzeug.utils import secure_filename

from redcap_preprocessing import redcap_preprocessing

app = Flask(__name__)

UPLOAD_FOLDER = os.path.join(os.path.dirname(__file__), 'uploads')
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

# Make sure the folder is empty
for filename in os.listdir(UPLOAD_FOLDER):
    file_path = os.path.join(UPLOAD_FOLDER, filename)
    try:
        if os.path.isfile(file_path):
            os.unlink(file_path)
    except Exception as e:
        print(e)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        # Check if a file part is present in the request
        if 'file' in request.files:
            file = request.files['file']
            if file.filename != '':
                filename = secure_filename(file.filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        
        # Fetch the selected analysis type
        analysis_type = request.form.get('analysis_type')

        redcap_preprocessing.preprocess_redcap_data(redcap_filepath=os.path.join(app.config['UPLOAD_FOLDER'], filename),
                             disease_type=analysis_type,
                             save_as_single_file=True)
                             
        return render_template('result.html')

    return render_template('index.html')

@app.route('/uploads/preprocessed_redcap_data/<filename>')
def download_file(filename):
    file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'preprocessed_redcap_data', filename)
    return send_file(file_path, as_attachment=True)

def run():
    app.run(debug=True)

if __name__ == '__main__':
    app.run(debug=True)
