# -*- coding:utf-8 -*-


from flask import Blueprint


AssClinical = Blueprint('AssClinical', __name__, template_folder="templates", static_folder='static')
from . import views
