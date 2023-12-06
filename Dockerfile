FROM python:3.8.5
COPY requirements.txt /python-flask/requirements.txt
ADD . /python-flask
WORKDIR /python-flask
RUN pip install -r requirements.txt

CMD [ "gunicorn", "-b 0.0.0.0:8000", "dashapp:server"]