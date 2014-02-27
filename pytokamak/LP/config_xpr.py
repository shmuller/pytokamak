from config import ShotNotFoundError

import config_xpr_200x
reload(config_xpr_200x)

import config_xpr_2012_spring
reload(config_xpr_2012_spring)

import config_xpr_2012_autumn
reload(config_xpr_2012_autumn)

import config_xpr_2013_spring
reload(config_xpr_2013_spring)

import config_xpr_2014_spring
reload(config_xpr_2014_spring)

campaign = config_xpr_200x.campaign \
         + config_xpr_2012_spring.campaign \
         + config_xpr_2012_autumn.campaign \
         + config_xpr_2013_spring.campaign \
         + config_xpr_2014_spring.campaign


