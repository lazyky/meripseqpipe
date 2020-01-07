from bottle import route, run
from bottle import get, post, request # or route
from bottle import error
from bottle import template
import os

@error(404)
def error404(error):
    return 'Nothing here, sorry'

@route('/')
@route('/hello/<name>')
def greet(name='Stranger'):
    return template('Hello {{name}}, how are you?', name=name)

@get('/login') # or @route('/login')
def login():
    return '''
        <form action="/login" method="post">
            Username: <input name="username" type="text" />
            Password: <input name="password" type="password" />
            <input value="Login" type="submit" />
        </form>
    '''

@post('/login') # or @route('/login', method='POST')
def do_login():
    username = request.forms.get('username')
    password = request.forms.get('password')
    if check_login(username, password):
        return "<p>Your login information was correct.</p>"
    else:
        return "<p>Login failed.</p>"

def check_login(username, password):
    if username == 'zky' and password == 'zky':
        return True
    else:
        return False



if __name__ == '__main__':
    results = os.system("ipconfig")
    print(results)

#run(debug=True, reloader=True)