# -*- coding:utf-8 -*-


from flask import Blueprint


zMap = Blueprint('zMap', __name__)
from . import views