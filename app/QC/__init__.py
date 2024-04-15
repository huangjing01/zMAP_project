# -*- coding:utf-8 -*-


from flask import Blueprint


QC = Blueprint('QC', __name__, template_folder="templates", static_folder='static')
from . import views
