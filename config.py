import os
import psycopg2

basedir = os.path.abspath(os.path.dirname(__file__))

def get_env_variable(name):
    try:
        return os.environ[name]
    except KeyError:
        message = "Expected environment variable '{}' not set.".format(name)
        raise Exception(message)




POSTGRES_URL="NA"
POSTGRES_USER="NA"
POSTGRES_PW="NA"
POSTGRES_DB="NA"

current_directory = os.getcwd()
DB_URL = 'sqlite:///'+ current_directory + '/your_database_name.db'
#B_URL = 'postgresql+psycopg2://{user}:{pw}@{url}/{db}'.format(user=POSTGRES_USER,pw=POSTGRES_PW,url=POSTGRES_URL,db=POSTGRES_DB)
#DB_URL = 'sqlite:////Users/tenzin/Desktop/projects/test/biodash/your_database_name.db'


#'DATABASE_URL'
class BaseConfig:
    SQLALCHEMY_DATABASE_URI = (DB_URL)
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    SECRET_KEY = 'secret_key'
    FLASK_ENV='development'
    #SECRET_KEY = os.environ['SECRET_KEY']


