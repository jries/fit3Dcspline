  msg = "The \"Cell Counter\" plugin has been replaced\n"
    + "by ImageJ's built in multi-point tool.";
  setTool("multipoint");
  showMessage("Cell Counter", msg);
  run("Point Tool...", "");
