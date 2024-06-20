#include <pthread.h>
#include <stdio.h>
#include <string>
#include <Python.h>
#include "MolTwisterRestAPI.h"

void CMolTwisterRestAPI::run() const
{
    pthread_t threadHandle;

    if(pthread_create(&threadHandle, nullptr, threadRun, nullptr))
    {
        printf("Error: Could not create REST API thread!\r\n");
    }
}

void* CMolTwisterRestAPI::threadRun(void*)
{
    std::string pythonString =
        R"(
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
    serve(app, host="127.0.0.1", port=5000)
    #app.run(port=5000) # Use this for debug output to console, instead of serve()
)";

    PyGILState_STATE state = PyGILState_Ensure();
    PyRun_SimpleString(pythonString.data());
    PyGILState_Release(state);
    return nullptr;
}
