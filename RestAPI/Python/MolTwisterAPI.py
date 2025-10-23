#
# Copyright (C) 2024 Richard Olsen.
# DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
#
# This file is part of MolTwister.
#
# MolTwister is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MolTwister is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
#

from flask import Flask, request
from flask_cors import CORS
from waitress import serve
import moltwister as mt

markdownHelpDoc = mt.get_help_as_markdown().replace('\r', '').split('\n')
tutorialCount = mt.get_tutorial_count()
tutorials = []

for i in range(0, tutorialCount):
    tutorials.append(mt.get_tutorial(i).replace('\r', '').split('\n'))

app = Flask(__name__)
CORS(app)

@app.get('/api_name')
def return_api_name():
    return {"api_name" : "MolTwister"}

@app.get('/startup_info')
def return_startup_info():
    return {"startup_info" : mt.get_startup_info()}

@app.get('/help')
def return_help():
    return {"help" : mt.get_help()}

@app.route('/help_cmd', methods=['POST'])
def return_help_cmd():
    jsonPayload = request.get_json(force=True)
    commandString = jsonPayload['cmd']
    commandLine = jsonPayload['cmdline']
    return {"help_cmd" : mt.get_help_for_command(commandString, commandLine)}

@app.get('/help_md')
def return_help_md():
    jsonPayload = request.get_json(force=True)
    if jsonPayload['ret'] == 'count':
        return {"help_md_cnt" : len(markdownHelpDoc)}
    if jsonPayload['ret'] == 'line':
        lineIndex = int(jsonPayload['index'])
        return {'help_md' : markdownHelpDoc[lineIndex]}

@app.get('/tutorial')
def return_tutorial():
    jsonPayload = request.get_json(force=True)
    if jsonPayload['ret'] == 'tutcount':
        return {"tutcnt" : len(tutorials)}
    if jsonPayload['ret'] == 'linecount':
        tutIndex = int(jsonPayload['tutindex'])
        return {"linecnt" : len(tutorials[tutIndex])}
    if jsonPayload['ret'] == 'line':
        tutIndex = int(jsonPayload['tutindex'])
        lineIndex = int(jsonPayload['lineindex'])
        return {'tutorial' : tutorials[tutIndex][lineIndex]}

@app.route('/exec', methods=['POST'])
def execute_command():
    jsonPayload = request.get_json(force=True)
    commandString = jsonPayload['cmd']
    return {"exec_ret" : mt.exec(commandString)}

@app.route('/query_atom_chunk_count', methods=['POST'])
def return_atom_chunk_count():
    jsonPayload = request.get_json(force=True)
    atomsInChunk = int(jsonPayload['atoms_in_chunk'])
    return {"atom_chunks" : mt.get_atom_chunk_count(atomsInChunk)}

@app.route('/query_atom_chunk', methods=['POST'])
def return_atom_chunk():
    jsonPayload = request.get_json(force=True)
    atomsInChunk = int(jsonPayload['atoms_in_chunk'])
    chunkIndex = int(jsonPayload['chunk_index'])
    return {"atom_chunk" : mt.get_atom_chunk(atomsInChunk, chunkIndex)}

@app.route('/reset_retrieved_chunks_state', methods=['POST'])
def reset_retrieved_chunks_state():
    return {"result" : "ok"}

@app.route('/query_deleted_atoms_chunk_count', methods=['POST'])
def return_deleted_atoms_chunk_count():
    jsonPayload = request.get_json(force=True)
    atomsInChunk = int(jsonPayload['atoms_in_chunk'])
    return {"atom_chunks" : mt.get_deleted_atoms_chunk_count(atomsInChunk)}

@app.route('/query_deleted_atoms_chunk_since_last_chunk_transfer', methods=['POST'])
def return_deleted_atoms_chunk_since_last_chunk_transfer():
    jsonPayload = request.get_json(force=True)
    atomsInChunk = int(jsonPayload['atoms_in_chunk'])
    chunkIndex = int(jsonPayload['chunk_index'])
    return {"deleted_chunk" : mt.get_deleted_atoms_chunk_since_last_chunk_transfer(atomsInChunk, chunkIndex)}

@app.route('/query_added_atoms_chunk_count', methods=['POST'])
def return_added_atoms_chunk_count():
    jsonPayload = request.get_json(force=True)
    atomsInChunk = int(jsonPayload['atoms_in_chunk'])
    return {"atom_chunks" : mt.get_added_atoms_chunk_count(atomsInChunk)}

@app.route('/query_added_atoms_chunk_since_last_chunk_transfer', methods=['POST'])
def return_added_atoms_chunk_since_last_chunk_transfer():
    jsonPayload = request.get_json(force=True)
    atomsInChunk = int(jsonPayload['atoms_in_chunk'])
    chunkIndex = int(jsonPayload['chunk_index'])
    return {"added_chunk" : mt.get_added_atoms_chunk_since_last_chunk_transfer(atomsInChunk, chunkIndex)}

if __name__ == "__main__":
    serve(app, host=mt.get_host_ip(), port=int(mt.get_host_port()))
    #app.run(port=5000) # Use this for debug output to console, instead of serve()
