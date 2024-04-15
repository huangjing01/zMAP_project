# -*- coding:utf-8 -*-


from flask import Blueprint


HVPannotation = Blueprint('HVPannotation', __name__, template_folder="templates", static_folder='static')
from . import views
