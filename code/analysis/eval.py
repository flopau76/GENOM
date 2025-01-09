import os
from typing import Tuple, List

from analysis.origin_search import Conclusion

def load_files(path_analysis_report : str, path_ref_report : str) -> Tuple[List[str], List[str]]:
    """
    Compare the analysis outputed file with the reference provided file.
    Verify only sending and receiving species, and the distance between 
    found distances with those predicted.
    """
    
    with open(path_analysis_report) as analysis:
        analysis = analysis.readlines()

    with open(path_ref_report) as report:
        report = report.readlines()
    
    report_list = [line.strip('\n').split('\t\t') for line in report]

    return analysis[1:], report_list[1:]


def compare_files(path_analysis_report : str, path_ref_report : str, window_size : int) -> Tuple[float, List[Conclusion]]:
    """
    Compare the analysis report with the reference of HGT occurences provided.
    """
    analysis, report = load_files(path_analysis_report, path_ref_report)

    successful_backtrack = 0
    valid = []
    for line in analysis:
        line = line.strip('\n')
        line_list = line.split('\t\t')

        liste_map = list(map(lambda x : x[0] == line_list[0] 
                                        and (int(line_list[1]) > int(x[1])-window_size and int(line_list[1]) < int(x[2])+window_size) 
                                        and (x[3] == line_list[2]),
                                        report))

        if True in liste_map:
            valid.append(
                Conclusion(
                    sender_found=line_list[0],
                    position_receiver=line_list[1],
                    receiver=line_list[2],
                    position_sender=line_list[3]
                    )
                )
            successful_backtrack+=1
    
    TP_rate = (successful_backtrack/len(report))*100
    return TP_rate, valid

"""    
print(compare_files(r'C:\Subpbiotech_cours\BT5\BIM_BMC\GENOM\project\project_git\GENOM\output\analysis\analysis_report.txt',
              r'C:\Subpbiotech_cours\BT5\BIM_BMC\GENOM\project\project_git\GENOM\output\output_generator\HGT_report.txt',
              5000))
"""