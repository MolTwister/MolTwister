from flask import Flask, request
from waitress import serve
import moltwister as mt

app = Flask(__name__)

@app.get('/api_name')
def return_api_name():
    return {"api_name" : "MolTwister"}

@app.route('/exec', methods=['POST'])
def execute_command():
    jsonPayload = request.get_json(force=True)
    commandString = jsonPayload['cmd']
    return mt.exec(commandString)

if __name__ == "__main__":
    serve(app, host="127.0.0.1", port=5000) # :TODO: Need to implement mt.get_host_ip() and mt.get_host_port()
    #app.run(port=5000) # Use this for debug output to console, instead of serve()
