import pytest
import sys
from app import create_app




@pytest.fixture
def app():
    flask_app = app.create_app()
    yield flask_app

    
@pytest.fixture
def client(app):
    return app.create_app()
