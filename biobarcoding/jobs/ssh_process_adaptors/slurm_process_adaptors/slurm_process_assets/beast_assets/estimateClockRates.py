import xml.etree.ElementTree as ET
import sys, re

def dict_input2dict(dict_input: str):
    result = re.findall(r"[^{}:,]+", dict_input)
    new_dict = {}
    for i in range(0, len(result), 2):
        new_dict[result[i].strip()] = result[i+1].strip()
    return new_dict

if __name__ == '__main__':
    LOCUS = dict_input2dict(sys.argv[1])
    print(LOCUS)
    method = sys.argv[2]
    MEANRATE_UPPER = {'chloroplast': "3.0E-9", 'ITS': "9.0E-9"}
    MEANRATE_LOWER = {'chloroplast': "1.0E-9", 'ITS': "1.0E-9"}
    MEANRATE_VALUE = {'chloroplast': "1.0E-9", 'ITS': "3.0E-9"}
    uniform_id = 1000

    UCLDMEAN_PARAM_TEMPLATE = '<parameter id="ucldMean.c:$id" spec="parameter.RealParameter" upper="$upper" lower="$lower" name="stateNode">$value</parameter>'
    MEANRATEPRIOR_PARAM_TEMPLATE = f"""     <prior id="MeanRatePrior.c:$id" name="distribution" x="@ucldMean.c:$id">
                        <Uniform id="Uniform.$uniform_id" name="distr" upper="$upper" lower="$lower" />
                    </prior>"""
    UCLMEANSCALER_PARAM_TEMPLATE = '<operator id="ucldMeanScaler.c:$id" spec="ScaleOperator" parameter="@ucldMean.c:$id" scaleFactor="0.5" weight="1.0"/>'
    RELAXEDUPDOWN_PARAM_TEMPLATE = '''<operator id="relaxedUpDownOperator.c:$id" spec="UpDownOperator" scaleFactor="0.75" weight="3.0">
                <up idref="ucldMean.c:$id"/>
                <down idref="Tree.t:treepartition"/>
            </operator>'''
    LOG_PARAM_TEMPLATE = '<log idref="ucldMean.c:$id"/>'

    root = ET.parse(f'beast_{method}_tmp.xml')

    partitions = []

    for data in root.iter('data'):
        if 'spec' in data.attrib.keys():
            partitions.append(data.attrib['id'])

    estimated_ids = []

    for branchRateModel in root.iter('branchRateModel'):
        tagId = branchRateModel.attrib['id'].split(':')[1]
        if not 'clock.rate' in branchRateModel.attrib.keys():
            branchRateModel.set('clock.rate', f'@ucldMean.c:{tagId}')
            for p in branchRateModel.iter('parameter'):
                if p.attrib['id'].startswith('ucldMean'):
                    branchRateModel.remove(p)
        else:
            estimated_ids.append(tagId)
            partitions.remove(tagId)

    for p in root.iter('parameter'):
        for estimated_id in estimated_ids:
            if 'id' in p.attrib.keys() and p.attrib['id'] == f'ucldMean.c:{estimated_id}':
                p.set('upper', MEANRATE_UPPER[LOCUS[estimated_id]])
                p.set('lower', MEANRATE_LOWER[LOCUS[estimated_id]])
                p.text = MEANRATE_VALUE[LOCUS[estimated_id]]

    for p in root.iter('prior'):
        for estimated_id in estimated_ids:
            if 'id' in p.attrib.keys() and p.attrib['id'] == f'MeanRatePrior.c:{estimated_id}':
                u = p.find('Uniform')
                u.set('upper', MEANRATE_UPPER[LOCUS[estimated_id]])
                u.set('lower', MEANRATE_LOWER[LOCUS[estimated_id]])

    for r in root.iter():
        if r.tag == 'state':
            state = r
            break

    for tagId in partitions:
        newParam = UCLDMEAN_PARAM_TEMPLATE.replace('$id', tagId).replace('$value', f'{MEANRATE_VALUE[LOCUS[tagId]]}').replace("$upper", f'{MEANRATE_UPPER[LOCUS[tagId]]}').replace("$lower", f'{MEANRATE_LOWER[LOCUS[tagId]]}')
        state.append(ET.fromstring(newParam))

    for distribution in root.iter('distribution'):
        if distribution.attrib['id'] == 'prior':
            for tagId in partitions:
                newPrior = MEANRATEPRIOR_PARAM_TEMPLATE.replace("$id", tagId).replace("$uniform_id", str(uniform_id)).replace("$upper", f'{MEANRATE_UPPER[LOCUS[tagId]]}').replace("$lower", f'{MEANRATE_LOWER[LOCUS[tagId]]}')
                distribution.append(ET.fromstring(newPrior))
                uniform_id += 1
    
    for r in root.iter():
        if r.tag == 'run':
            run = r
            break

    for operator in run.iterfind('operator'):
        if 'id' in operator.attrib.keys() and operator.attrib['id'] == 'FixMeanMutationRatesOperator':
            run.remove(operator)

    for tagId in partitions:
        newOperator = UCLMEANSCALER_PARAM_TEMPLATE.replace("$id", tagId)
        run.append(ET.fromstring(newOperator))
        newOperator = RELAXEDUPDOWN_PARAM_TEMPLATE.replace("$id", tagId)
        run.append(ET.fromstring(newOperator))

    for logger in root.iter('logger'):
        if logger.attrib['id'] == 'tracelog':
            for tagId in partitions:
                newLog = LOG_PARAM_TEMPLATE.replace("$id", tagId)
                logger.append(ET.fromstring(newLog))
            break

    for r in root.iter():
        if r.tag == 'beast':
            r.set('beautistatus', 'noAutoSetClockRate')
            break

    root.write(f"beast_{method}.xml")
