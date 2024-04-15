# -*- coding:utf-8 -*-


from flask import Blueprint


Sample_Subgroup = Blueprint('Sample_Subgroup', __name__)
from . import views
