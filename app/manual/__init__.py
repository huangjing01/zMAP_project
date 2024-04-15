# -*- coding:utf-8 -*-


from flask import Blueprint


manual = Blueprint('manual', __name__)
from . import views