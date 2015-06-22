Search.setIndex({envversion:47,filenames:["authors","contributing","eventbuilder","faq","format","history","index","modules","pax","pax.plugins","pax.plugins.corrections","pax.plugins.for_tests","pax.plugins.io","pax.plugins.posrec","pax.plugins.signal_processing","plugin","position_reconstruction","pyroot","readme","usage"],objects:{"":{pax:[8,0,0,"-"]},"pax.FolderIO":{InputFromFolder:[8,3,1,""],WriteToFolder:[8,3,1,""]},"pax.FolderIO.InputFromFolder":{close:[8,2,1,""],file_extension:[8,4,1,""],get_all_events_in_current_file:[8,2,1,""],get_events:[8,2,1,""],get_first_and_last_event_number:[8,2,1,""],get_single_event:[8,2,1,""],get_single_event_in_current_file:[8,2,1,""],init_file:[8,2,1,""],open:[8,2,1,""],select_file:[8,2,1,""],shutdown:[8,2,1,""],startup:[8,2,1,""]},"pax.FolderIO.WriteToFolder":{close:[8,2,1,""],close_current_file:[8,2,1,""],file_extension:[8,4,1,""],open:[8,2,1,""],open_new_file:[8,2,1,""],shutdown:[8,2,1,""],startup:[8,2,1,""],write_event:[8,2,1,""],write_event_to_current_file:[8,2,1,""]},"pax.core":{Processor:[8,3,1,""]},"pax.core.Processor":{fallback_configuration:[8,4,1,""],get_metadata:[8,2,1,""],get_plugin_by_name:[8,2,1,""],get_plugin_search_paths:[8,7,1,""],instantiate_plugin:[8,2,1,""],load_configuration:[8,2,1,""],process_event:[8,2,1,""],run:[8,2,1,""],setup_logging:[8,2,1,""],shutdown:[8,2,1,""]},"pax.data_model":{ListField:[8,3,1,""],Model:[8,3,1,""],StrictModel:[8,3,1,""]},"pax.data_model.Model":{from_bson:[8,1,1,""],from_json:[8,1,1,""],get_fields_data:[8,2,1,""],get_list_field_info:[8,1,1,""],to_bson:[8,2,1,""],to_dict:[8,2,1,""],to_json:[8,2,1,""]},"pax.datastructure":{Event:[8,3,1,""],Hit:[8,3,1,""],Peak:[8,3,1,""],Pulse:[8,3,1,""],ReconstructedPosition:[8,3,1,""],SumWaveform:[8,3,1,""]},"pax.datastructure.Event":{S1s:[8,2,1,""],S2s:[8,2,1,""],all_hits:[8,4,1,""],dataset_name:[8,4,1,""],duration:[8,2,1,""],empty_event:[8,1,1,""],event_number:[8,4,1,""],get_sum_waveform:[8,2,1,""],get_sum_waveform_names:[8,2,1,""],is_channel_suspicious:[8,4,1,""],length:[8,2,1,""],n_channels:[8,4,1,""],n_hits_rejected:[8,4,1,""],noise_pulses_in:[8,4,1,""],peaks:[8,4,1,""],pulses:[8,4,1,""],sample_duration:[8,4,1,""],start_time:[8,4,1,""],stop_time:[8,4,1,""],sum_waveforms:[8,4,1,""]},"pax.datastructure.Hit":{area:[8,4,1,""],center:[8,4,1,""],channel:[8,4,1,""],found_in_pulse:[8,4,1,""],height:[8,4,1,""],index_of_maximum:[8,4,1,""],is_rejected:[8,4,1,""],left:[8,4,1,""],length:[8,4,1,""],noise_sigma:[8,4,1,""],right:[8,4,1,""]},"pax.datastructure.Peak":{area:[8,4,1,""],area_fraction_top:[8,4,1,""],area_per_channel:[8,4,1,""],bottom_hitpattern_spread:[8,4,1,""],contributing_channels:[8,4,1,""],detector:[8,4,1,""],does_channel_contribute:[8,4,1,""],does_channel_have_noise:[8,4,1,""],height:[8,4,1,""],hit_time_mean:[8,4,1,""],hit_time_std:[8,4,1,""],hits:[8,4,1,""],index_of_maximum:[8,4,1,""],left:[8,4,1,""],mean_amplitude_to_noise:[8,4,1,""],n_contributing_channels:[8,4,1,""],n_noise_channels:[8,4,1,""],noise_channels:[8,4,1,""],range_50p_area:[8,4,1,""],range_90p_area:[8,4,1,""],reconstructed_positions:[8,4,1,""],right:[8,4,1,""],top_hitpattern_spread:[8,4,1,""],type:[8,4,1,""]},"pax.datastructure.Pulse":{baseline:[8,4,1,""],channel:[8,4,1,""],left:[8,4,1,""],length:[8,4,1,""],maximum:[8,4,1,""],minimum:[8,4,1,""],noise_sigma:[8,4,1,""],raw_data:[8,4,1,""],right:[8,4,1,""]},"pax.datastructure.ReconstructedPosition":{algorithm:[8,4,1,""],goodness_of_fit:[8,4,1,""],ndf:[8,4,1,""],phi:[8,4,1,""],r:[8,4,1,""],x:[8,4,1,""],y:[8,4,1,""]},"pax.datastructure.SumWaveform":{channel_list:[8,4,1,""],detector:[8,4,1,""],is_filtered:[8,2,1,""],name:[8,4,1,""],name_of_filter:[8,4,1,""],samples:[8,4,1,""]},"pax.event_builder":{run:[8,5,1,""]},"pax.exceptions":{OutputFileAlreadyExistsError:[8,6,1,""],PulseBeyondEventError:[8,6,1,""]},"pax.formats":{BulkOutputFormat:[8,3,1,""],HDF5Dump:[8,3,1,""],NumpyDump:[8,3,1,""],PandasCSV:[8,3,1,""],PandasFormat:[8,3,1,""],PandasHDF5:[8,3,1,""],PandasHTML:[8,3,1,""],PandasJSON:[8,3,1,""],ROOTDump:[8,3,1,""]},"pax.formats.BulkOutputFormat":{close:[8,2,1,""],data_types_present:[8,4,1,""],file_extension:[8,4,1,""],open:[8,2,1,""],prefers_python_strings:[8,4,1,""],read_data:[8,2,1,""],supports_append:[8,4,1,""],supports_array_fields:[8,4,1,""],supports_read_back:[8,4,1,""],supports_write_in_chunks:[8,4,1,""],write_data:[8,2,1,""]},"pax.formats.HDF5Dump":{close:[8,2,1,""],data_types_present:[8,4,1,""],file_extension:[8,4,1,""],n_in_data:[8,2,1,""],open:[8,2,1,""],read_data:[8,2,1,""],supports_array_fields:[8,4,1,""],supports_read_back:[8,4,1,""],supports_write_in_chunks:[8,4,1,""],write_data:[8,2,1,""]},"pax.formats.NumpyDump":{close:[8,2,1,""],data_types_present:[8,4,1,""],f:[8,4,1,""],file_extension:[8,4,1,""],n_in_data:[8,2,1,""],open:[8,2,1,""],read_data:[8,2,1,""],supports_array_fields:[8,4,1,""],supports_read_back:[8,4,1,""],write_data:[8,2,1,""]},"pax.formats.PandasCSV":{pandas_format_key:[8,4,1,""]},"pax.formats.PandasFormat":{open:[8,2,1,""],pandas_format_key:[8,4,1,""],supports_array_fields:[8,4,1,""],write_data:[8,2,1,""],write_pandas_dataframe:[8,2,1,""]},"pax.formats.PandasHDF5":{close:[8,2,1,""],data_types_present:[8,4,1,""],file_extension:[8,4,1,""],n_in_data:[8,2,1,""],open:[8,2,1,""],prefers_python_strings:[8,4,1,""],read_data:[8,2,1,""],string_data_length:[8,4,1,""],supports_append:[8,4,1,""],supports_read_back:[8,4,1,""],supports_write_in_chunks:[8,4,1,""],write_pandas_dataframe:[8,2,1,""]},"pax.formats.PandasHTML":{pandas_format_key:[8,4,1,""]},"pax.formats.PandasJSON":{pandas_format_key:[8,4,1,""]},"pax.formats.ROOTDump":{close:[8,2,1,""],data_types_present:[8,4,1,""],file_extension:[8,4,1,""],n_in_data:[8,2,1,""],numpy_type:[8,4,1,""],open:[8,2,1,""],read_data:[8,2,1,""],root_type:[8,4,1,""],supports_array_fields:[8,4,1,""],supports_read_back:[8,4,1,""],supports_write_in_chunks:[8,4,1,""],write_data:[8,2,1,""]},"pax.plugin":{BasePlugin:[8,3,1,""],InputPlugin:[8,3,1,""],OutputPlugin:[8,3,1,""],TransformPlugin:[8,3,1,""]},"pax.plugin.BasePlugin":{process_event:[8,2,1,""],shutdown:[8,2,1,""],startup:[8,2,1,""]},"pax.plugin.InputPlugin":{get_events:[8,2,1,""],get_single_event:[8,2,1,""],number_of_events:[8,4,1,""],process_event:[8,2,1,""]},"pax.plugin.OutputPlugin":{process_event:[8,2,1,""],write_event:[8,2,1,""]},"pax.plugin.TransformPlugin":{process_event:[8,2,1,""],transform_event:[8,2,1,""]},"pax.plugins":{corrections:[10,0,0,"-"],for_tests:[11,0,0,"-"],io:[12,0,0,"-"],posrec:[13,0,0,"-"],signal_processing:[14,0,0,"-"]},"pax.plugins.corrections":{ExampleCorrection:[10,0,0,"-"]},"pax.plugins.corrections.ExampleCorrection":{ExampleCorrection:[10,3,1,""]},"pax.plugins.corrections.ExampleCorrection.ExampleCorrection":{startup:[10,2,1,""],transform_event:[10,2,1,""]},"pax.plugins.for_tests":{Dummy:[11,0,0,"-"]},"pax.plugins.for_tests.Dummy":{DummyInput:[11,3,1,""],DummyOutput:[11,3,1,""],DummyTransform2:[11,3,1,""],DummyTransform3:[11,3,1,""],DummyTransform:[11,3,1,""]},"pax.plugins.for_tests.Dummy.DummyInput":{get_events:[11,2,1,""]},"pax.plugins.for_tests.Dummy.DummyOutput":{write_event:[11,2,1,""]},"pax.plugins.io":{Avro:[12,0,0,"-"],BSON:[12,0,0,"-"],BulkOutput:[12,0,0,"-"],MongoDB:[12,0,0,"-"],Pickle:[12,0,0,"-"],RawWaveformDump:[12,0,0,"-"],WaveformSimulator:[12,0,0,"-"],XED:[12,0,0,"-"]},"pax.plugins.io.Avro":{ReadAvro:[12,3,1,""],WriteAvro:[12,3,1,""]},"pax.plugins.io.Avro.ReadAvro":{close:[12,2,1,""],file_extension:[12,4,1,""],get_all_events_in_current_file:[12,2,1,""],open:[12,2,1,""]},"pax.plugins.io.Avro.WriteAvro":{close:[12,2,1,""],file_extension:[12,4,1,""],open:[12,2,1,""],startup:[12,2,1,""],write_event_to_current_file:[12,2,1,""]},"pax.plugins.io.BSON":{ReadBSON:[12,3,1,""],ReadJSON:[12,3,1,""],ReadZippedBSON:[12,3,1,""],WriteBSON:[12,3,1,""],WriteJSON:[12,3,1,""],WriteZippedBSON:[12,3,1,""]},"pax.plugins.io.BSON.ReadBSON":{close:[12,2,1,""],file_extension:[12,4,1,""],get_all_events_in_current_file:[12,2,1,""],open:[12,2,1,""]},"pax.plugins.io.BSON.ReadJSON":{close:[12,2,1,""],file_extension:[12,4,1,""],get_all_events_in_current_file:[12,2,1,""],open:[12,2,1,""]},"pax.plugins.io.BSON.ReadZippedBSON":{close:[12,2,1,""],file_extension:[12,4,1,""],get_single_event_in_current_file:[12,2,1,""],open:[12,2,1,""]},"pax.plugins.io.BSON.WriteBSON":{close:[12,2,1,""],file_extension:[12,4,1,""],open:[12,2,1,""],write_event_to_current_file:[12,2,1,""]},"pax.plugins.io.BSON.WriteJSON":{close:[12,2,1,""],file_extension:[12,4,1,""],open:[12,2,1,""],write_event_to_current_file:[12,2,1,""]},"pax.plugins.io.BSON.WriteZippedBSON":{close:[12,2,1,""],file_extension:[12,4,1,""],open:[12,2,1,""],write_event_to_current_file:[12,2,1,""]},"pax.plugins.io.BulkOutput":{BulkOutput:[12,3,1,""],ReadFromBulkOutput:[12,3,1,""]},"pax.plugins.io.BulkOutput.BulkOutput":{get_index_of:[12,2,1,""],shutdown:[12,2,1,""],startup:[12,2,1,""],write_event:[12,2,1,""]},"pax.plugins.io.BulkOutput.ReadFromBulkOutput":{convert_record:[12,2,1,""],get_events:[12,2,1,""],startup:[12,2,1,""]},"pax.plugins.io.MongoDB":{IOMongoDB:[12,3,1,""],MongoDBReadUntriggered:[12,3,1,""],MongoDBReadUntriggeredFiller:[12,3,1,""],MongoDBWriteTriggered:[12,3,1,""],sampletime_fmt:[12,5,1,""]},"pax.plugins.io.MongoDB.IOMongoDB":{chunks:[12,7,1,""],number_events:[12,2,1,""],setup_access:[12,2,1,""],setup_input:[12,2,1,""],shutdown:[12,2,1,""],startup:[12,2,1,""],update_run_doc:[12,2,1,""]},"pax.plugins.io.MongoDB.MongoDBReadUntriggered":{get_events:[12,2,1,""],sliding_window2:[12,2,1,""],sliding_window:[12,2,1,""],startup:[12,2,1,""]},"pax.plugins.io.MongoDB.MongoDBReadUntriggeredFiller":{startup:[12,2,1,""],transform_event:[12,2,1,""]},"pax.plugins.io.MongoDB.MongoDBWriteTriggered":{startup:[12,2,1,""],write_event:[12,2,1,""]},"pax.plugins.io.Pickle":{DirWithPickleFiles:[12,3,1,""],ReadFromStackedPickleFolder:[12,3,1,""],WriteToPickleFile:[12,3,1,""],WriteToStackedPickleFolder:[12,3,1,""]},"pax.plugins.io.Pickle.DirWithPickleFiles":{get_events:[12,2,1,""],get_single_event:[12,2,1,""],startup:[12,2,1,""]},"pax.plugins.io.Pickle.ReadFromStackedPickleFolder":{close:[12,2,1,""],file_extension:[12,4,1,""],get_all_events_in_current_file:[12,2,1,""],open:[12,2,1,""]},"pax.plugins.io.Pickle.WriteToPickleFile":{write_event:[12,2,1,""]},"pax.plugins.io.Pickle.WriteToStackedPickleFolder":{close:[12,2,1,""],file_extension:[12,4,1,""],open:[12,2,1,""],write_event_to_current_file:[12,2,1,""]},"pax.plugins.io.RawWaveformDump":{DumpSumWaveformToBinary:[12,3,1,""],WaveformDumperBase:[12,3,1,""]},"pax.plugins.io.RawWaveformDump.DumpSumWaveformToBinary":{get_waveform_to_dump:[12,2,1,""]},"pax.plugins.io.RawWaveformDump.WaveformDumperBase":{get_waveform_to_dump:[12,2,1,""],startup:[12,2,1,""],write_event:[12,2,1,""]},"pax.plugins.io.WaveformSimulator":{WaveformSimulator:[12,3,1,""],WaveformSimulatorFromCSV:[12,3,1,""],WaveformSimulatorFromNEST:[12,3,1,""]},"pax.plugins.io.WaveformSimulator.WaveformSimulator":{get_events:[12,2,1,""],get_instructions_for_next_event:[12,2,1,""],s1:[12,2,1,""],s2:[12,2,1,""],shutdown:[12,2,1,""],simulate_single_event:[12,2,1,""],startup:[12,2,1,""],store_true_peak:[12,2,1,""]},"pax.plugins.io.WaveformSimulator.WaveformSimulatorFromCSV":{get_instructions_for_next_event:[12,2,1,""],shutdown:[12,2,1,""],startup:[12,2,1,""]},"pax.plugins.io.WaveformSimulator.WaveformSimulatorFromNEST":{get_instructions_for_next_event:[12,2,1,""],startup:[12,2,1,""],variables:[12,4,1,""]},"pax.plugins.io.XED":{XedInput:[12,3,1,""]},"pax.plugins.io.XED.XedInput":{close:[12,2,1,""],file_extension:[12,4,1,""],get_first_and_last_event_number:[12,2,1,""],get_single_event_in_current_file:[12,2,1,""],open:[12,2,1,""]},"pax.plugins.posrec":{NeuralNet:[13,0,0,"-"],PosRecChiSquareGamma:[13,0,0,"-"],PosSimple:[13,0,0,"-"]},"pax.plugins.posrec.NeuralNet":{NeuralNet:[13,3,1,""],PosRecNeuralNet:[13,3,1,""]},"pax.plugins.posrec.NeuralNet.NeuralNet":{run:[13,2,1,""],run_layer:[13,2,1,""]},"pax.plugins.posrec.NeuralNet.PosRecNeuralNet":{startup:[13,2,1,""],transform_event:[13,2,1,""]},"pax.plugins.posrec.PosRecChiSquareGamma":{PosRecChiSquareGamma:[13,3,1,""]},"pax.plugins.posrec.PosRecChiSquareGamma.PosRecChiSquareGamma":{function_chi_square_gamma:[13,2,1,""],shutdown:[13,2,1,""],startup:[13,2,1,""],transform_event:[13,2,1,""]},"pax.plugins.posrec.PosSimple":{PosRecWeightedSum:[13,3,1,""]},"pax.plugins.posrec.PosSimple.PosRecWeightedSum":{startup:[13,2,1,""],transform_event:[13,2,1,""]},"pax.plugins.signal_processing":{CheckPulses:[14,0,0,"-"],Cluster:[14,0,0,"-"],HitFinder:[14,0,0,"-"]},"pax.plugins.signal_processing.CheckPulses":{CheckBounds:[14,3,1,""],ConcatenateAdjacentPulses:[14,3,1,""],SortPulses:[14,3,1,""]},"pax.plugins.signal_processing.CheckPulses.CheckBounds":{startup:[14,2,1,""],transform_event:[14,2,1,""]},"pax.plugins.signal_processing.CheckPulses.ConcatenateAdjacentPulses":{transform_event:[14,2,1,""]},"pax.plugins.signal_processing.CheckPulses.SortPulses":{transform_event:[14,2,1,""]},"pax.plugins.signal_processing.Cluster":{ClusterPlugin:[14,3,1,""],GapSize:[14,3,1,""],HitDifference:[14,3,1,""],MeanShift:[14,3,1,""]},"pax.plugins.signal_processing.Cluster.ClusterPlugin":{cluster_hits:[14,2,1,""],startup:[14,2,1,""],transform_event:[14,2,1,""]},"pax.plugins.signal_processing.Cluster.GapSize":{cluster_hits:[14,2,1,""],startup:[14,2,1,""]},"pax.plugins.signal_processing.Cluster.HitDifference":{cluster_hits:[14,2,1,""]},"pax.plugins.signal_processing.Cluster.MeanShift":{cluster_hits:[14,2,1,""],get_gap_probability:[14,7,1,""],get_gap_size:[14,2,1,""],startup:[14,2,1,""]},"pax.plugins.signal_processing.HitFinder":{FindHits:[14,3,1,""],compute_hit_properties:[14,4,1,""],compute_pulse_properties:[14,4,1,""],find_intervals_above_threshold:[14,4,1,""]},"pax.plugins.signal_processing.HitFinder.FindHits":{startup:[14,2,1,""],transform_event:[14,2,1,""]},"pax.simulation":{SimulatedHitpattern:[8,3,1,""],Simulator:[8,3,1,""],truncated_gauss_rvs:[8,5,1,""]},"pax.simulation.Simulator":{get_luminescence_times:[8,2,1,""],make_hitpattern:[8,2,1,""],pmt_pulse_current:[8,2,1,""],s1_photons:[8,2,1,""],s2_electrons:[8,2,1,""],s2_scintillation:[8,2,1,""],singlet_triplet_delays:[8,2,1,""],to_pax_event:[8,2,1,""]},"pax.utils":{InterpolatingMap:[8,3,1,""],Memoize:[8,3,1,""],chunk_in_ntuples:[8,5,1,""],cluster_by_diff:[8,5,1,""],data_file_name:[8,5,1,""],get_named_configuration_options:[8,5,1,""],mad:[8,5,1,""],weighted_mean_variance:[8,5,1,""]},"pax.utils.InterpolatingMap":{data_field_names:[8,4,1,""],get_value:[8,2,1,""],get_value_at:[8,2,1,""],plot:[8,2,1,""]},pax:{FolderIO:[8,0,0,"-"],core:[8,0,0,"-"],data_model:[8,0,0,"-"],datastructure:[8,0,0,"-"],event_builder:[8,0,0,"-"],exceptions:[8,0,0,"-"],formats:[8,0,0,"-"],plugin:[8,0,0,"-"],plugins:[9,0,0,"-"],simulation:[8,0,0,"-"],units:[8,0,0,"-"],utils:[8,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","classmethod","Python class method"],"2":["py","method","Python method"],"3":["py","class","Python class"],"4":["py","attribute","Python attribute"],"5":["py","function","Python function"],"6":["py","exception","Python exception"],"7":["py","staticmethod","Python static method"]},objtypes:{"0":"py:module","1":"py:classmethod","2":"py:method","3":"py:class","4":"py:attribute","5":"py:function","6":"py:exception","7":"py:staticmethod"},terms:{"2to3":17,"3d_pl":19,"3d_plotter":19,"__1":17,"__dict__":17,"__import__":17,"__pycache__":8,"__version__":3,"_base":19,"_importhook":17,"_intern":8,"break":1,"case":[8,12,15,17],"class":[4,5,6,8,10,11,12,13,14],"default":[4,5,8,12,15,18],"export":[3,18],"final":12,"float":[4,8],"function":[1,5,8,12,13,14,15,16,17],"import":[3,5,8,17],"int":[4,8,14],"long":12,"new":[1,4,5,8,12,14,16],"return":[4,8,12,13,14,15,17],"static":[5,8,12,14],"super":18,"switch":[5,18],"true":[3,4,8,12,16,18],"try":[3,8,12,14,17,18],"while":15,a_buss09:0,aalber:0,abbrevi:8,abl:[3,8,18],about:[1,3,4,8,12,14],abov:[3,8,14,17,18],absolut:8,access:[8,12,18],account:14,accuraci:16,action:8,activ:[13,17],actual:1,adc:[4,5,8,14],add:[1,3,8],addit:[8,15],addition:13,address:[3,12],adjac:14,advanc:19,advis:15,after:[3,5,12,14,15,17,19],again:3,agreement:5,alex:13,algorithm:[4,5,8,13,14,16],all:[4,5,8,12,13,15],all_hit:[4,8],allow:[5,14,18],almost:8,alon:16,along:5,alper:0,alreadi:[1,14,16],also:[1,4,5,8,12,15,16,17,18],altern:[17,18],alwai:[1,3,4,5,8,12,15],among:17,amplitud:[4,8],anaconda3:[3,18],anaconda:[1,3,18],analys:17,analyz:4,ani:[1,4,8,12,14,15,18],annot:18,anoth:[1,16],another_map_nam:8,anybodi:1,anyth:[1,12,14],apach:12,apache_avro:12,apart:8,api:12,append:[8,12,13,16],append_data:12,applic:[3,8],appreci:1,appropri:12,apt:17,aquisit:12,arbitrari:5,arcan:8,archiv:18,area:[4,8,12,14],area_fraction_top:[4,8],area_per_channel:[4,8],aren:1,arg:[8,12],argmax:14,argument:[1,4,5,8,18,19],around:[4,5],arrai:[4,8,12,13,14,16],arriv:[8,12],ascend:[4,8],ask:1,assign:1,associ:4,assum:[12,14,17,18],assumpt:16,atanh:13,attribut:8,author:1,automat:[1,3,4,14,15,18],avail:[12,15,19],averag:13,avro:[5,8,9],awar:4,axel:0,axi:8,backbon:8,backend:[5,8,12],backward:17,bad:5,bandpass:5,bart:0,bartp:[0,3],base:[5,8,10,11,12,13,14,15,16,18],baselin:[8,14],baselineexcursionmethod:5,baseplugin:8,bash:18,bashrc:3,basic:[5,14],becaus:[3,4,8,12],becom:8,been:19,befor:1,behaviour:[5,17],belong:[4,8],below:[4,8,14],bern:[5,19],best:1,better:[3,5,17],between:[1,4,5,8,12,14],beyond:14,bia:13,bias:13,bin:[3,4,8,15,17,18,19],binari:5,bind:5,binutil:17,bit:[1,4,8],blah:8,blahblah:14,block:5,bool:[4,8],both:16,bottom:[4,5,8,16,18],bottom_hitpattern_spread:[4,8],bound:[4,8,14],boundari:[8,14],branch:[1,8],breanch:1,breur:0,browser:18,bsd:18,bson:[5,8,9],buffer_s:12,bugfix:[1,5],build:[1,5,8,12,17],builtin:17,bulk:[5,8],bulkoutput:[8,9],bulkoutputformat:[8,12],bunch:12,buß:0,bypass:17,bzip2:12,calcul:[13,16],call:[4,8,12,16,17],can:1,cannot:[3,8,18],capabl:5,care:8,cartesian:8,caus:[8,17],center:[4,8,14],centermost:[4,8],central:4,centroid:13,cern:17,certif:3,cflag:[3,18],chain:[5,12,16],chang:[1,5,8,12,17,19],channel:[4,5,8,14],channel_list:[4,8],channel_waveform:12,channels_to_use_for_reconstruct:16,charact:15,charg:[5,13],check:[1,3,4,8,12,14,16,17],checkbound:14,checkout:1,checkpuls:[8,9],chi:[5,13],chi_square_gamma:13,child:12,chose:1,christoph:[0,3],chunk:[8,12],chunk_in_ntupl:8,chunk_siz:8,class_nam:15,class_to_load_to:12,classif:[5,12],classmethod:[4,8],clean:5,clean_shutdown:8,cleanli:15,cleanup:[5,17],clear:1,clearli:15,click:3,clock:[4,8],clone:[1,3],close:[8,12,15],close_current_fil:8,cluster:[4,5,8,9,12],cluster_by_diff:8,cluster_hit:14,clusterplugin:14,code1:8,code:[1,4,5,12,15,17,18],code_qu:1,coderr:0,coincid:[4,8],collect:12,column:12,com:[1,3,8,18],come:[1,12],comic:1,command:[3,5,17,18],comment:[5,14,15,17,19],commit:[1,5],common:12,commonli:[4,8],commun:4,comparison:5,compat:17,compil:[3,17,18],complain:3,complet:5,compress:[12,18],comput:[4,5,8,12,14],compute_hit_properti:14,compute_pulse_properti:14,concaten:[12,14],concatenateadjacentpuls:14,conda:[3,18],confid:1,config:[4,8,14,15,17,19],config_dict:8,config_nam:8,config_path:[8,18,19],config_str:8,config_to_init:8,config_valu:[8,10,11,12,13,14],configur:[4,5,8,12,15,16,17,18,19],connect:13,consid:1,consist:8,constant:15,construct:16,constructor:[4,8,15],contact:1,contain:[4,5,8,12,14,15,16,17],content:[6,7],continu:5,continuum:18,contributing_channel:[4,8],control:5,convert:[4,5,8,12,14,17],convert_numpy_arrays_to:8,convert_record:12,convolut:5,coordin:8,coordinate_system:8,copi:[1,3,15,17],core:[5,7],correct:[5,8,9],correctli:17,correspond:12,cosin:[5,18],could:[1,4,8,17,18],count:[4,8,14,16],cout:17,cover:14,coverag:5,coveral:5,creat:[1,4,5,6,8],creation:8,csv:8,ctunnel:0,culpa:5,currenli:16,current:[5,8,12,16],custom:15,cut:5,cylindr:8,cython:18,dai:1,daniel:0,daq:[5,8,12,18],daq_injector:19,daqinjector:5,darkmatt:13,data:[3,4,5,8,12,14,15,16,17,18,19],data_field_nam:8,data_file_nam:8,data_model:[4,7],data_types_pres:8,databas:[5,12,15,19],datafram:5,dataset:[5,12],dataset_nam:[4,8],datastructur:[4,5,7],datatyp:8,date:1,deactiv:17,dead:5,deal:[16,17],debug:[5,14,19],decai:8,decid:14,declar:8,decor:14,def:17,defin:[5,8,15],definit:[5,8,17],degre:16,delet:8,demo:5,depend:[3,4,5,8,14,15,17,18],deposit:8,deprec:[4,8,12],depth:[8,12],deriv:15,descend:[4,8],describ:[8,15,16],descript:[1,8,15,19],deseri:[4,8],desktop:17,destruct:5,destructor:15,detail:[1,12,13],detector:[4,5,8],determin:[5,8,12,13],dev:[17,18],deviat:[4,8],df_name:8,diagnost:8,dict:8,dictionari:[5,8,15],diff_threshold:8,differ:[5,14,15],difficult:[1,12],digit:[4,5,8,12,18],dir:8,directli:[1,14],directori:[5,8,12,15,16,17,18,19],dirwithpicklefil:12,disabl:17,discuss:1,disk:[12,18],displai:18,distribut:[13,16,18],divis:8,dll:17,dname:12,doc:[1,5,8],docstr:[1,8,12],doe:[1,4,8,12,14,15],does_channel_contribut:[4,8],does_channel_have_nois:[4,8],don:[1,8,14],done:[1,3,4,8,16],down:[5,8,17],download:17,dpkg:17,drift:[8,12],drupal:17,dsp:5,dtype:[4,8],due:[5,12],dummi:[8,9],dummyinput:11,dummyoutput:11,dummytransform2:11,dummytransform3:11,dummytransform:11,dump:[8,12],dumpsumwaveformtobinari:12,durat:[4,8],each:[4,5,8,12,13,15,16],earli:1,earliest:[4,8],eas:[4,8,18],easier:[1,18],easili:5,easy_instal:3,echo:17,edg:[4,5,8],edit:8,effici:5,either:[12,14,16,17,18],electron:[8,12],electron_arrival_tim:8,electron_tim:12,electrons_gener:8,element_typ:8,elr:8,els:[8,12],email:[3,18],empti:[4,5,8],empty_ev:[4,8],emul:[5,8],enabl:[5,14,17],encod:[5,12],end:[5,8,12,17],energi:8,enforc:8,enhanc:5,ensur:[1,8,14],enter:3,entir:[5,12],entiti:12,env:17,environ:[1,3,17],epoch:8,error:[12,14,17],etc:[8,12,15],evalu:8,even:[1,4,8],evenli:8,event_build:7,event_class:12,event_numb:[4,8],event_posit:[8,12],eventu:12,ever:12,everi:[1,4,8,12,15,17],everyth:[3,8],exactli:3,exampl:[4,5,8,12,16,17,18],examplecorrect:[8,9],exce:14,except:[5,7],excim:8,execut:3,exim:8,exist:[1,8,12,13,16],exit:[17,19],expand:[1,4],expect:3,experi:[3,12,18],explain:[1,4,8],explor:[1,18],extend:[5,8,14],extens:[1,5,12],extern:5,extra:[5,12],extra_path:8,fact:3,fail:3,fallback_configur:8,fals:[4,8,18],faq:18,fast:12,faster:5,favorit:8,favourit:3,fax:[5,8,12],feed:13,feel:[1,3,4],few:8,fiction:3,field:[5,8,12,15],fields_to_ignor:[8,12],file:[1,3,5,8,12,15,16,17,18,19],file_extens:[8,12],filenam:[8,12],fill:[8,14],fillvalu:8,filter:[4,5,8,18],find:[1,5,8,13,14,16,18],find_intervals_above_threshold:14,findbigpeak:5,finder:[5,12],findhit:14,fine:12,first:[4,5,6,8,12,17],first_ev:8,first_event_numb:8,fit:[4,8,13,16],flag:[1,4,8],flake8:1,flat:12,float32:8,float64:[4,8],folder:[1,5,8,12],folderio:7,follow:[1,3,4,8,12,15,16,17,18],food:8,for_test:[8,9],forget:14,fork:[1,5],form:16,format:3,forward:13,found:[4,8,12,14,15,17],found_in_puls:[4,8],four:17,foward:13,fraction:[4,8],framework:[5,8,12],free:[1,3,4,18],freedom:16,frequenc:5,frequent:1,from:[3,4,5,8,12,13,14,15,17,18,19],from_bson:8,from_json:8,full:16,fulli:12,fun:5,function_chi_square_gamma:13,further:5,gain:[8,13],gamma:[5,13],gap:14,gap_width:14,gapsiz:[5,14],gate:8,gauss:8,gcc:[3,17],gener:[1,4,5,8,12,16],generate_mock_correction_map:8,get_all_events_in_current_fil:[8,12],get_data_field:8,get_ev:[8,11,12,15],get_fields_data:8,get_first_and_last_event_numb:[8,12],get_gap_prob:14,get_gap_s:14,get_index_of:12,get_instructions_for_next_ev:12,get_list_field_info:8,get_luminescence_tim:8,get_metadata:8,get_named_configuration_opt:8,get_plugin_by_nam:8,get_plugin_search_path:8,get_single_ev:[8,12],get_single_event_in_current_fil:[8,12],get_sum_waveform:[4,8],get_sum_waveform_nam:[4,8],get_valu:8,get_value_at:8,get_waveform_to_dump:12,getattr:17,git:1,git_ssl_no_verifi:3,github:[1,3,13,18],give:12,given:[4,8],glob:8,global:16,gnu:17,goat:3,gohlk:3,good:[4,8,13,16,17],goodness_of_fit:[4,8],googl:1,googlecod:18,graviti:[4,8,14],greatli:1,grid:8,group:[5,13,14],groupbaudi:13,guess:3,guid:[1,15],guillaum:12,gzip:[12,17],h5py:18,hadoop:12,handi:18,handl:[5,8],happen:1,happi:17,hard:15,hasattr:17,hater:3,have:[3,4,5,8,12,15,17,18,19],hdf5:[5,8],hdf5dump:8,head:8,height:[4,8],help:[1,5,18,19],helper:8,henc:12,here:[1,4,8,12,18],hidden:13,high:[4,5,8],high_threshold:14,higher:[5,12],highest:[4,8],hint:18,hit_time_mean:[4,8],hit_time_std:[4,8],hitdiffer:14,hitfind:[5,8,9],hitlist:8,hitpattern:[4,8],hoc:17,hold:[4,8],home:[3,17,18],hour:1,how:1,howev:[1,12,18],html:[1,8,13],http:[1,3,8,12,13,14,17,18],idem:8,ident:[5,13],ignor:3,img:1,important_modul:3,imprecis:[4,8],improv:5,impuls:5,inc:17,incdir:17,includ:[1,3,4,5,8,15,16,17,18],inclus:[4,8,14],increment:1,index:[4,5,6,8,12,14],index_of_maximum:[4,8],indic:1,individu:[4,8,12,14],info:18,inform:[4,8,12,13,16,18],inherit:5,ini:[5,15,16,18,19],init:8,init_fil:8,initi:[8,13,15],initial_baseline_sampl:14,initializi:8,inject:5,injector:5,input:[5,8,12,13],input_valu:13,inputfromfold:[8,12],inputplugin:[8,11,12,15],insist:3,instal:[1,3,5,6,8],instanc:[4,8,12],instantiate_plugin:8,instruct:[3,5,12,18],int16:8,int32:8,int64:[4,8],integ:[4,8],integr:[4,5,8,12],intend:12,interact:[8,12,16,17,18],interfac:[5,8,12],intermedi:15,intern:12,interpol:[5,8],interpolatingmap:8,interpret:5,interv:14,intput:15,intro:6,invers:14,iomongodb:12,ioniz:[4,8],is_channel_suspici:[4,8],is_filt:[4,8],is_reject:[4,8],isn:[8,18],isol:5,issu:[1,3,5],iter:[1,8,12],ith:8,jaalber:0,januari:[4,8],jell:0,jit:14,json:[8,12],just:[1,3,5,18],just_test:8,kaminski:3,keep:1,kei:[4,8],kept:12,kind:12,kish:13,kish_thesiselectron:13,know:[3,18],kodiaq:12,kwarg:[4,8],kwargs_dict:[4,8],label:14,laptop:17,larg:[1,12],largest:5,last:[4,8,11,12,14,17],last_ev:[8,11],lastli:[1,4],later:[4,5,8,14],latest:17,launch:5,layer:[13,18],lce:[13,16],lce_map_fil:16,lce_map_file_nam:16,ldflag:[3,18],leak:5,learn:18,least:[4,8],leav:16,left:[4,8,12,14],left_boundari:8,leftmostleft:[4,8],length:[4,5,8,12,14],let:[3,18],level:[4,5,8,12,14,19],lexic:8,lhep:0,lib:[3,17,18],libdir:17,libpyroot:17,libpython3:17,libpython:17,librari:[8,12,17,18],libsnappi:18,libx11:17,libxdio:12,libxext:17,libxft:17,libxpm:17,licens:18,lifetim:8,like:[3,5,15,18],line:[3,5,17,18],linearli:8,lint:1,linux:[17,18],linux_setup:3,list:[1,4,8,12,15,16,17,18],listfield:[4,8],littl:1,load:[5,8,13,18],load_configur:8,local:[1,3],locat:[8,15,17,18],log:[5,8,19],log_scale_entire_ev:18,logger:[5,8],login:3,lone_puls:14,longer:[12,14],look:[1,3,17,18],lost:1,lot:5,low:16,low_threshold:14,lowest:13,macport:3,mad:8,made:8,mai:[1,3,8,12],maintain:[1,8],major:[1,5,8,17],make:[1,3,4,8,15,17,18,19],make_hitpattern:8,manag:[1,18],mani:[1,5],manipul:8,manual:17,map:[5,8,13,16],map_nam:8,mask:5,master:1,match:[4,5,8],matplotlib:18,matrix:8,matur:5,max:[5,14],max_differ:14,max_gap_s:14,maxima:14,maximum:[4,8,12,14],mea:5,mean:[4,5,8,12,14,16],mean_amplitude_to_nois:[4,8],mean_shift:14,meaning:8,meaningless:5,meanshift:[5,14],meant:8,meantim:3,median:8,meet:1,member:15,memoiz:8,memori:[5,12],mention:3,merg:[1,5],mesh:8,mess:8,messag:19,metadata:12,method:[5,8,14,15,16],micromodel:5,middl:[12,14],might:1,min:14,min_width:5,minim:13,minimum:[8,13,16],minor:[5,8,17],minut:1,miss:[3,5],mkdir:17,mname:12,mode:[3,5,8,16],model:8,modif:[15,17],modifi:[5,8,14,15,16,18],modified_root_v5:17,modified_root_v6:17,modul:[3,5,6,7],moment:12,mongo:[5,19],mongodb:[5,8,9],mongodbreaduntrigg:12,mongodbreaduntriggeredfil:12,mongodbwritetrigg:12,more:[1,3,5,8,12,14,18,19],most:[3,5,15,16,18],motiv:5,move:5,much:[3,5],muenster:0,multi:17,multipl:[5,12],music:5,must:[4,5,8,12,15,18],mutil:5,my_fil:18,my_mean:8,my_postprocess:15,my_std:8,mynewclassthatdoessomethingtransform:15,n_channel:[4,8],n_contributing_channel:[4,8],n_hidden:13,n_hits_reject:[4,8],n_in_data:8,n_input:13,n_noise_channel:[4,8],n_output:13,n_photon:8,n_rv:8,n_x:8,nag:3,name:[1,4,8,12,15,16,19],name_of_filt:[4,8],nan:8,nanosecond:[4,8],narrow:1,nativ:14,natur:14,navig:18,ndarrai:8,ndf:[4,8,13,16],nearli:[5,15],need:[1,3,4,8,12,15,16,17],neg:[5,8],nest:[5,8],nest_i:12,nest_nel:12,nest_nph:12,nest_nr:12,nest_t:12,nest_x:12,nest_z:12,net:[5,13,18],network:13,neural:[5,13,18],neuralnet:[8,9],neuron:13,newdsp:[5,19],newlin:12,next:[14,18],nice:[3,8,18],nikhef:0,no_reconstruct:16,nobodi:1,nois:[4,5,8,14],noise_channel:[4,8],noise_pulses_in:[4,8],noise_sigma:[4,8,14],none:[4,8,12],nonneg:8,nonstandard:18,note:[4,8,17],noth:14,now:[1,3,5,17,18],npz:8,nth:12,nuclear:12,num:12,numba:[5,14,18],number:[4,8,12,14,16,17],number_ev:12,number_of_ev:8,numer:[5,15],numpi:[8,12,13,14,18],numpy_typ:8,numpydump:8,object:[4,8,12,13,15,16],occur:[5,8,12],off:[4,5,8],offer:3,offici:1,offset:[8,12],often:[1,18],old:[3,5,8,12],older:18,onc:[1,8,12,18],onli:[3,4,8,12,15,16],onlin:[5,18],only_reconstruct:16,open:[1,5,8,12,15],open_new_fil:8,oper:1,opportun:5,opt:3,optim:[12,16],option:[6,12],order:[8,15,16],ordinari:[4,8],org:[8,12,14],origin:[1,8,16],other:[1,3,5,8,12,13,15,16,18],otherwis:[4,8,14],otter:19,our:[5,8,12,18],out:[5,8,12,17],output:[3,4,5,8,12,13],output_dir:18,output_format:12,output_nam:12,outputfilealreadyexistserror:8,outputplugin:[8,11,12,15],over:[8,12],overhang:5,overhead:12,overload:15,overrid:[5,8,15],overwrit:12,overwrite_data:12,own:[15,18,19],p_value_offset:14,packag:7,pad:8,page:[3,6],panda:[5,18],pandas_format_kei:8,pandascsv:8,pandasformat:8,pandashdf5:8,pandashtml:8,pandasjson:8,parallel:[5,8],param:[8,14],paramet:[4,8,12,13,16],parent_configur:18,pars:5,part:[1,17],partial:[4,8],particl:[4,8],particular:[12,19],pass:[1,4,8,13],past:5,patch:8,path:[3,8,15,18,19],pattern:[5,13,16,18],pax:1,pax_dir:8,pax_info:12,paxer:[3,18,19],paxit:[5,15],pdf:13,peak_typ:12,peakfind:5,pelsser:0,peopl:[3,8,18],pep8:1,per:[4,8,12],perform:[5,15],permiss:18,person:1,phd:12,phi:[4,8],photoelectron:[4,8],photon:[5,8,12],photon_tim:[8,12],physic:[4,8,12],physik:13,pickl:[5,8,9],pickler:5,pip:[1,3,17,18],place:8,plai:18,plane:16,plant:12,pleas:[1,3,4,8,15,18],plot:[3,5,8,18,19],plot_to_dir:19,plotchannelwaveforms3d:18,ploteventsummari:18,plug:12,plugin:[5,6,7],plugin_path:[8,15],pmt:[4,8,12,13,16,18],pmt_pulse_curr:8,png:1,point:[3,4,8,12,17,18],poisson:16,polish:5,port:[5,12],posit:[4,5,8,13,15],posrec:[8,9],posrecchi:16,posrecchisquaregamma:[4,8,9],posrecneuralnet:13,posrecweightedsum:[13,16],possibl:[1,12,14,16],possimpl:[5,8,9],powel:16,practis:17,pre:17,prefers_python_str:8,prefix:18,prerequisit:17,present:12,previou:[5,13],print:[3,17],privat:18,probabl:[14,15],problem:3,proce:1,procedur:12,process:[5,8,12,15,16,17,18,19],process_ev:8,processor:5,produc:[12,18],product:8,program:[17,18],project:[1,5,12],promis:8,proper:[8,15,17],properti:[5,8,12,14],propos:1,provid:[4,8,12,16,17],puls:[4,5,8,12,14],pulsebeyondeventerror:8,pure:14,purpos:[4,8],push:1,put:[1,8,17],py34:17,py3:[5,17],pymongo3:5,pyroot:3,pytabl:18,python2:8,python3:[3,8,17],python:[1,3,4,5,8,14,17,18],pythondir:17,quick:8,quickli:14,quit:14,rais:[5,12,14,18],random:8,rang:[4,8,12],range_50p_area:[4,8],range_90p_area:[4,8],ranges_buffer_s:12,rather:[4,8,17],ratio:8,raw:[5,8,12,14,18],raw_data:8,raw_data_fil:8,raw_hit:14,rawwaveformdump:[8,9],reach:8,reactiv:18,read:[5,8,12,19],read_data:8,readavro:12,readbson:12,readfrombulkoutput:12,readfromstackedpicklefold:12,readi:[1,17],readjson:12,readm:1,readzippedbson:12,realli:[1,8],reason:[4,8],reboot:3,recalcul:14,receiv:12,recoil:12,recoil_typ:[8,12],recommend:[1,18],reconstruct:[4,5,8,13],reconstructed_posit:[4,8],record:[4,8,12],recurs:8,redefin:15,redo:12,refactor:5,refer:8,region:[5,8,14],regular:8,reject:[4,5,8],rel:[8,14],relat:[3,5],releas:[4,5,8],relev:18,remot:12,remov:5,renam:8,reorgan:5,repeat:[1,5],repect:15,replac:12,repo:18,repositori:18,reprocess:[5,12],reproduc:1,requir:[3,5,6,15,17],resolv:1,respect:[8,17],respons:[5,8,12],rest:12,restart:12,restrict:8,result:[1,3,4,8,12,14],result_buff:14,return_indic:8,revers:[4,8],review:1,right:[4,8,12,14],right_boundari:8,rightmost:[4,8],rightmostright:[4,8],root_:17,root_typ:8,rootdump:8,round:[1,4,8],routin:[8,12],rst:1,run:1,run_14:18,run_lay:13,runtim:12,s1_default_recombination_tim:12,s1_photon:[8,12],s2_electron:[8,12],s2_scintil:8,s2_size:14,s2_width:14,sacrific:3,sad:3,safe:8,sai:8,said:3,same:[4,8,12,14],sampl:[4,8,12,14],sample_dur:[4,8],sampletime_fmt:12,sander:0,sanderb:0,save:[12,15,19],scalar:8,scale:12,scheme:12,scientif:18,scikit:18,scintil:[4,8],scipi:[8,16,18],scope:1,screen:19,script:[3,5,15],search:[6,8,12,14,15,17],second:8,secondli:13,section:[8,15,18],secur:3,see:[3,4,5,8,12,13,16,18,19],seen:5,select:8,select_fil:8,self:[8,11,15,17],send:1,separ:[5,12,14,18],sequenti:12,serial:[4,8,12],server:5,set:1,setup:[1,3,5,18],setup_access:12,setup_input:12,setup_log:8,sever:[3,5,8,12,18],shape:8,shard:5,shell:3,shift:14,ship:[4,8],should:[1,3,4,5,8,14,15,17,18],shouldn:3,show:[3,4,8,17,19],shut:[5,8,17],shutdown:[8,12,13,15],side:5,sigma:[4,8],signal:[4,5,8,13,16,18],signal_process:[8,9],signific:[4,8],simpl:13,simpledsp:5,simplifi:5,simul:[5,7],simulate_single_ev:12,simulatedhitpattern:8,sinc:[4,5,8,16,17],singl:[5,8,12,13],singlet:8,singlet_ratio:8,singlet_triplet_delai:8,site:8,size:[8,12],skeleton:8,slide:12,sliding_window2:12,sliding_window:12,slow:[12,14],small:[1,5,8,12],smaller:1,softwar:[3,5,18],solut:17,solv:[5,18],some:[3,5,8,12,14,17,18],some_field:8,somebodi:1,someth:[1,4,8,18],sometyp:8,somewher:[12,14],sophist:18,sort:[4,8,14],sort_kei:[4,8],sortpuls:14,sourc:[3,4,5,8,10,11,12,13,14,15,17,18],space:8,spe:14,special:17,specif:[1,6,8],specifi:[4,5,8,12,15,17],speed:[5,12,16],spend:1,sphinx:1,spit:8,split:8,squar:[4,5,8,13],src:17,stabl:[5,8],stackedpickl:12,stackoverflow:8,stai:8,standard:[8,17],start_tim:[4,8,12],startup:[8,10,12,13,14,15],stat:8,state:8,statement:17,statist:[5,16],std:[4,8,14,17],step:[1,3,6,15,17],still:[1,8,14],stolen:8,stop:[4,5,8,14,19],stop_aft:19,stop_tim:[4,8],stoptim:8,storag:12,store:[4,8,11,12],store_true_peak:12,straightforward:8,strang:5,stream:[4,8],strictmodel:8,string:12,string_data_length:[8,12],structer:8,structur:[5,8,12],stuck:3,student:0,stuff:[8,12],style:1,subclass:8,subdir:8,subdirectori:15,submodul:7,subpackag:7,subplot:5,subset:1,success:[12,17],sudo:17,sum:[4,5,8,12,13],sum_waveform:[4,8],summari:5,sumwaveform:[4,8],support:[5,8,12],supports_append:8,supports_array_field:8,supports_read_back:8,supports_write_in_chunk:8,sure:[3,15,17],suspici:[4,5,8],syntax:8,system:[1,3,8,12,18],tag:1,take:[1,8,12,14,16],taken:16,talk:8,tar:[17,18],tbranch:8,team:3,temporari:[4,8],test:[1,3,4,5,8,11,18],test_pax:1,than:[1,3,5,8,12,14,18],thei:[1,5,8,17,18],them:[1,3,14,15,17],therefor:[1,3,4,8,16],thesi:[12,13],thi:[1,3,4,5,8,12,13,14,15,16,17,18,19],thing:[5,18],think:[4,8],thisroot:17,those:[3,12,13],though:[12,18],thought:8,thread:17,threshold:[4,8],through:[1,3,8,18],thu:16,time:[4,8,12,14,15,17],timestamp:[8,12],to_bson:8,to_dict:8,to_fil:8,to_json:8,to_pax_ev:8,todo:12,tofil:12,too:[4,8,12],top:[3,4,5,8,13,16,18],top_hitpattern_spread:[4,8],topcuoglu:0,topçuoğlu:0,total:[8,12],tox:1,tpc:[4,8],tradit:5,transform:[5,8],transform_ev:[8,10,12,13,14,15],transformplugin:[8,10,11,12,13,14,15],travi:[1,3,5,17],tree:[4,8],tri:8,trick:8,trigger:[4,5,8,12],triplet:8,troubl:8,troubleshoot:1,truncat:[8,12,14],truncated_gauss_rv:8,truncnorm:8,truth:[5,12],ttree:8,tunnel:[0,3],tupl:8,tweak:17,two:[3,5,15,16],txt:5,typic:[1,18],ubuntu:[3,17,18],uint16:[4,8],uncom:18,uncompress:12,underestim:16,understand:1,undon:1,unexpectedli:1,unfortun:3,uni:0,unib:0,unit:[4,7],unittest:1,unix:[4,8],unknown:[4,8,14],unless:[4,8],unpack:17,untest:1,until:5,untrigg:12,updat:[1,3,4,5,18],update_run_doc:12,updatedb:17,upon:4,usag:6,user:[8,15,18],usr:17,usual:8,util:7,uva:0,uzh:13,v1724:[4,8],valu:[5,8,12,13,15,16],valuex1y1:8,valuex1y2:8,valuex2y1:8,valuex2y2:8,variabl:[8,12,15,17],variou:[4,5,8],veri:[1,12,13],version:[1,3,5,12,14,17],veto:[4,5,8,18],via:[1,15,16,18],view:5,virtual:[1,17],virtualenv:[1,17],virualenv:17,wai:[1,4,8],want:[1,8,12,14,15,17,18,19],warn:3,wave:8,waveform:3,waveformdumperbas:12,waveformsimul:[5,8,9],waveformsimulatorfromcsv:12,waveformsimulatorfromnest:12,web:[5,18],webpag:18,websit:[1,3],weight:[4,5,8,13],weighted_mean_vari:8,welcom:1,well:[12,16],were:5,wget:18,what:[1,3,4,8,12,15,18],whatev:1,when:[1,3,5,8,12,17],whenev:[8,12],where:[4,8,12,17],whether:1,which:[3,4,5,8,12,13,14,15,16,18],white:5,whl:3,who:[1,3,4,8],whoever:1,whose:[4,8],width:5,wiki:[1,12,14],wikipedia:[12,14],win:3,within:[3,4,8,12,15,18],without:[4,8,12,13],won:14,work:1,workaround:3,would:[1,8,12,19],wrap:5,write_data:8,write_ev:[8,11,12,15],write_event_to_current_fil:[8,12],write_in_chunk:12,write_pandas_datafram:8,writeavro:12,writebson:12,writejson:12,writepanda:5,writetofold:[8,12],writetopicklefil:12,writetostackedpicklefold:12,writezippedbson:12,written:[1,12],www:13,x86_64:[17,18],x_max:8,x_min:8,xam:[5,19],xdio:12,xdp:[4,8],xe100_150213_1411:18,xe100_150213_1411_000000:18,xed:[5,8,9],xedinput:[12,18],xenon100:[5,8,12,18,19],xenon1:19,xenon1t:[1,3,5],xenon:[3,4,8,13,17,18],xerawdp:[4,5,8],xkcd:1,xml:5,xvfz:18,yield:[12,13,15,16],you:[1,3,4,8,12,14,15,17,18,19],your:[1,3,8,15,18,19],your_map_nam:8,yourself:[3,17],zero:[5,8],zip:12,zipfil:12,zle0:12,zle:5},titles:["Credits","Contributing","Event builder","Frequently Asked Questions","Event format","History","Welcome to processor for Analyzing XENON1T&#8217;s documentation!","pax","pax package","pax.plugins package","pax.plugins.corrections package","pax.plugins.for_tests package","pax.plugins.io package","pax.plugins.posrec package","pax.plugins.signal_processing package","Plugins","Position Reconstruction","PyROOT installation","Installation","Usage"],titleterms:{"class":15,analyz:6,ask:3,avro:12,bson:12,bug:1,builder:2,bulkoutput:12,can:3,charg:16,checkpuls:14,chi:16,cluster:14,content:[8,9,10,11,12,13,14],contribut:1,contributor:0,core:8,correct:10,creat:15,credit:0,data_model:8,datastructur:8,develop:[0,1],document:[1,6],dummi:11,event:[2,4],event_build:8,examplecorrect:10,except:8,featur:[1,18],feedback:1,first:18,fix:1,folderio:[8,12],for_test:11,format:[4,8],frequent:3,gamma:16,get:[1,3],git:3,guidelin:1,histori:5,hit:4,hitfind:14,how:3,implement:1,indic:6,input:15,instal:[17,18],intro:[15,16],lead:0,lng:3,machin:3,minim:16,modul:[8,9,10,11,12,13,14],mongodb:12,neuralnet:13,nikhef:3,option:15,osx:3,output:15,packag:[8,9,10,11,12,13,14],pax:[3,7,8,9,10,11,12,13,14,18],peak:4,pickl:12,plugin:[8,9,10,11,12,13,14,15],posit:16,posrec:13,posrecchisquaregamma:13,possimpl:13,processor:6,pull:1,pyroot:17,question:3,rawwaveformdump:12,reconstruct:16,reconstructedposit:4,report:1,request:1,requir:18,root:[3,12],run:3,set:[3,16],signal_process:14,simul:8,snappi:3,specif:15,squar:16,start:1,step:18,submit:1,submodul:[8,10,11,12,13,14],subpackag:[8,9],sum:16,tabl:6,tip:1,transform:15,type:1,unit:8,usag:[16,19],util:8,waveform:4,waveformsimul:12,weight:16,welcom:6,window:3,work:3,write:1,xeclust:3,xed:12,xenon1t:6}})