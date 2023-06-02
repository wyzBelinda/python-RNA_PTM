import os.path


def read_cfg(cfg_path):
    import configparser

    cfg = configparser.ConfigParser()
    cfg.read(cfg_path)

    return cfg


def get_cfg_items(cfg, section_name):
    # 获取所有的配置项
    items = cfg.items(section_name)
    return dict(items)


def backup_cfgs(out_dir, cfg):
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, 'cfg.ini'), 'w') as configfile:
        cfg.write(configfile)

    import shutil
    shutil.copy(cfg['paths']['element_path'], os.path.join(out_dir, 'element.ini'))
    shutil.copy(cfg['paths']['icon_path'], os.path.join(out_dir, 'nts.ini'))
    shutil.copy(cfg['paths']['enzyme_path'], os.path.join(out_dir, 'enzyme.csv'))
    shutil.copy(cfg['paths']['mods_path'], os.path.join(out_dir, 'mods.csv'))
    shutil.copy(cfg['paths']['known_mod_path'], os.path.join(out_dir, 'considered_mods.txt'))


if __name__ == '__main__':
    config = read_cfg('data/cfg.ini')

    # items_dict = get_cfg_items(config, 'digest')
    # print(items_dict)

    # # 读取cfg文件
    # config = read_cfg('data/cfg.ini')
    #
    # paths_dict = get_cfg_items(config, "paths")

    # 备份参数文件
    backup_cfgs(config['paths']['out_dir'],config)
