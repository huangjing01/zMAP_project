# -*- coding:utf-8 -*-


from flask import Blueprint


AssMutation = Blueprint('AssMutation', __name__, template_folder="templates", static_folder='static')
from . import views
