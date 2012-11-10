import config_xpr_2012_spring
reload(config_xpr_2012_spring)

import config_xpr_2012_autumn 
reload(config_xpr_2012_autumn)

campaign = config_xpr_2012_spring.campaign \
         + config_xpr_2012_autumn.campaign


